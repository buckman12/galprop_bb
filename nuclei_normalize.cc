
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nuclei_normalize.cc *                         galprop package * 08/16/2001
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// renormalization coefficient is defined from PRIMARY protons flux only;
// this renormalization coefficient then applied to all nuclei and secondary species

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"
#include"ErrorLogger.h"
#include <cstring>
#include <sstream>

int Galprop::nuclei_normalize()
{
	INFO("Entry");
	ostringstream buf;

// identify the primary CR protons                            IMOS20000609
	int iprotons=-1;
	for(int i=n_species-1; i>=0; i--) if(gcr[i].Z==1&&gcr[i].A==1) {
			iprotons=i;      // IMOS20010816
			break;
		}
	if(iprotons==-1) {
		WARNING("CR protons not found!");
		return 1;
	}
	if(strcmp(gcr[iprotons].name,"Hydrogen_1")!=0) {
		WARNING("CR protons not found!");
		return 1;
	}
	buf.str("");
	buf<<"  CR protons found as species #"<<iprotons;
	INFO(buf.str());

	double v1,v2,v3,v4,v5,v6;
	double r0=0.5; // solar Galactocentric radius, kpc

	//   int ip=(int)(log(galdef.proton_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 0.5);//IMOS20060420
	//   int iz=(int)((-galdef.z_min)/galdef.dz + 0.5);//IMOS20060420

	// corrected here and elsewhere: want grid points containg normalizing point for interpolation
	int ip=(int)(log(galdef.proton_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 1.0e-6);//AWS20131107
	int iz=(int)((-galdef.z_min)                             /galdef.dz               + 1.0e-6);//AWS20131107


// check normalization point in Ekin  AWS20131107 NB galaxy.n_pgrid not defined!
	if(ip<0 || ip+1>gcr[iprotons].n_pgrid-1) {
		buf.str("");
		buf<<"ERROR: nuclei normalization Ekin: ip="<<ip<<": "<< galdef.proton_norm_Ekin << " outside range of grid "
		   <<  gcr[iprotons].Ekin[0] <<" - "<<  gcr[iprotons].Ekin[gcr[iprotons].n_pgrid-1]
		   <<" gcr[iprotons].n_pgrid="<<gcr[iprotons].n_pgrid
		   << " ; no normalization can be done, exiting!";
		INFO(buf.str());
		exit(1);
	}



	if(galdef.n_spatial_dimensions==2) {
		int ir=(int)((r0-galdef.r_min)/galdef.dr + 1.0e-6);//AWS20131107
		buf.str("");
		buf<<"Grid point for normalization: ir r[ir] iz z[iz] ip Ekin[ip]:Ekin[ip+1 "<<ir<<" " <<gcr[iprotons].r[ir]
		   <<" " <<iz <<" "<<gcr[iprotons].z[iz]
		   <<" "<<ip<<" "<< gcr[iprotons].Ekin[ip]<<":"<< gcr[iprotons].Ekin[ip+1]
		   <<" normalization Ekin="<<galdef.proton_norm_Ekin ;   //AWS20131107
		INFO(buf.str());

		v1=gcr[iprotons].cr_density.d2[ir  ][iz].s[ip];
		v2=gcr[iprotons].cr_density.d2[ir+1][iz].s[ip];
		v3=gcr[iprotons].cr_density.d2[ir  ][iz].s[ip+1];
		v4=gcr[iprotons].cr_density.d2[ir+1][iz].s[ip+1];
		v5=v1+(r0-gcr[iprotons].r[ir])/galdef.dr*(v2-v1); // r0 ip
		v6=v3+(r0-gcr[iprotons].r[ir])/galdef.dr*(v4-v3); // r0 ip+1
	}

	if(galdef.n_spatial_dimensions==3) {
		int ix=(int)((r0-galdef.x_min)/galdef.dx + 1.0e-6);//AWS20131107
		int iy=(int)(( 0-galdef.y_min)/galdef.dy + 1.0e-6);//AWS20131107

		buf.str("");
		buf<<"Grid point for normalization: ix x[ix] iy y[iy] iz z[iz] ip Ekin[ip]:Ekin[ip+1] "<<ix<<" " <<gcr[iprotons].x[ix]
		   <<" "<<iy<<" "<< gcr[iprotons].y[iy]<<" "<<iz <<" "<<gcr[iprotons].z[iz]
		   <<" "<<ip<<" "<< gcr[iprotons].Ekin[ip]<<":"<< gcr[iprotons].Ekin[ip+1]
		   <<" normalization Ekin="<<galdef.proton_norm_Ekin ;   //AWS20131107
		INFO(buf.str());

		v1=gcr[iprotons].cr_density.d3[ix  ][iy][iz].s[ip];     //AWS20001121
		v2=gcr[iprotons].cr_density.d3[ix+1][iy][iz].s[ip];     //AWS20001121
		v3=gcr[iprotons].cr_density.d3[ix  ][iy][iz].s[ip+1];   //AWS20001121
		v4=gcr[iprotons].cr_density.d3[ix+1][iy][iz].s[ip+1];   //AWS20001121
		v5=v1+(r0-gcr[iprotons].x[ix])/galdef.dx*(v2-v1); // r0 ip
		v6=v3+(r0-gcr[iprotons].x[ix])/galdef.dx*(v4-v3); // r0 ip+1
	}

	double vnorm=exp( log(v5)+log(galdef.proton_norm_Ekin/gcr[iprotons].Ekin[ip])/log(galdef.Ekin_factor)*log(v6/v5) );
	buf.str("");
	buf<<"v1 v2 v3 v4 v5 v6 vnorm  "<<v1<<" " <<v2<<" " <<v3 <<" "<<v4 <<" "<<v5<<" "<< v6<<" "<<vnorm;
	INFO(buf.str());

	galdef.source_normalization *= galdef.proton_norm_flux/vnorm; // IMOS20030214
//   cout<<" nuclei_normalize >>> source_normalization= "<<galdef.source_normalization<<endl;

// normalize all species except primary electrons and primary positrons since these are normalized independently
	for(int i=0; i<n_species; i++) {
		if(strcmp(gcr[i].name,"primary_electrons")==0) continue;// IMOS20050912
		if(strcmp(gcr[i].name,"primary_positrons")==0) continue;//  AWS20131023

// don't re-normalize the DM annihilation products
		if(strcmp(gcr[i].name,"DM_positrons"  )==0) continue;  // IMOS20050912
		if(strcmp(gcr[i].name,"DM_electrons"  )==0) continue;  // IMOS20050912
		if(strcmp(gcr[i].name,"DM_antiprotons")==0) continue;  // IMOS20050912
		if(strcmp(gcr[i].name,"DM_antideuterons")==0) continue;  // IMOS20050912
		if(strcmp(gcr[i].name,"DM_antihelium")==0) continue;  // IMOS20050912

		buf.str("");
		buf<<" normalizing "<<gcr[i].name;
		INFO(buf.str());
		if(galdef.verbose>=1) {
			buf.str("");
			buf<<" nucleus "<<gcr[i].name<<" before normalization:";
			INFO(buf.str());
			gcr[i].cr_density.print();
		}

		gcr[i].cr_density         *=(galdef.proton_norm_flux/vnorm);
		gcr[i].normalization_factor=(galdef.proton_norm_flux/vnorm);//AWS20010121

		if(galdef.verbose>=1) {
			buf.str("");
			buf<<" nucleus "<<gcr[i].name<<" after normalization:";
			INFO(buf.str());
			gcr[i].cr_density.print();
		}

	}
	if(galdef.verbose>=1) {
		INFO("primary protons after normalization:");
		gcr[iprotons].cr_density.print();
	}
	INFO("Exit");
	return 0;
}
