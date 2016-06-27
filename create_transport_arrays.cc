
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_transport_arrays.cc *                  galprop package * 10/12/2003
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>
#include "ErrorLogger.h"

#include <string>
#include <sstream>

static int key=-1;

int Galprop::create_transport_arrays(Particle &particle)
{
	INFO("Entry");
	int stat=0, A1,Z2,A2,K_electron;                                               // IMOS20010816
	int galdef_network_par=0;         // imos network, use everywhere                 IMOS20010816
	double t_half[2];                                                              // IMOS20010816
	double fragment_p,fragment_He,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann; // IMOS20000607

// ASSIGNING PRIMARY SOURCE FUNCTION

	INFO("Assigning primary source function");

	double spec_shape;

	double
	g_0      =galdef.nuc_g_0,             //AWS20131111
	g_0_inner=galdef.nuc_g_0,
	rigid_br0=galdef.nuc_rigid_br0,       //AWS20131111
	g_1      =galdef.nuc_g_1,
	g_1_inner=galdef.nuc_g_1_inner,
	rigid_br =galdef.nuc_rigid_br,
	g_2      =galdef.nuc_g_2,
	g_2_inner=galdef.nuc_g_2_inner,
	rigid_br2=galdef.nuc_rigid_br2,       //AWS20131111
	g_3      =galdef.nuc_g_3;             //AWS20131111



	const string priElecStr = "primary_electrons";

	if (priElecStr == particle.name) {
		g_0      =galdef.electron_g_0;                                    // IMOS20031012
		g_0_inner=galdef.electron_g_0_inner;                              // IMOS20031012
		rigid_br0=galdef.electron_rigid_br0;
		g_1_inner=galdef.electron_g_1_inner;
		g_1      =galdef.electron_g_1;
		rigid_br =galdef.electron_rigid_br;
		g_2      =galdef.electron_g_2;
		rigid_br2=galdef.electron_rigid_br2;  //AWS20131111
		g_3      =galdef.electron_g_3;        //AWS20131111
	}

	const string priPosiStr = "primary_positrons";                      //AWS20131023

	if (priPosiStr == particle.name) {                                  //AWS20131023
		g_0      =galdef.positron_g_0;
		rigid_br0=galdef.positron_rigid_br0;
		g_1      =galdef.positron_g_1;
		rigid_br =galdef.positron_rigid_br;
		g_2=      galdef.positron_g_2;
		rigid_br2=galdef.positron_rigid_br2;  //AWS20131111
		g_3      =galdef.positron_g_3;        //AWS20131111
	}

	if (g_0_inner == -100)
		g_0_inner = g_0;
	if (g_1_inner == -100)
		g_1_inner = g_1;
	if (g_2_inner == -100)
		g_2_inner = g_2;


	ostringstream ost;
	ost<<particle.name<<" g_0="<<g_0<<"  rigid_br0= "<<rigid_br0  // IMOS20031012
	   <<" g_1="<<g_1<<"  rigid_br= " <<rigid_br <<" g_2="<<g_2
	   <<" g_0_inner="<<g_0_inner << " g_1_inner="<<g_1_inner;
	INFO(ost.str());

	particle.primary_source_function = 0.; // IMOS20020418 whole 2D/3D particle.primary_source_function loops are changed

	//To allow for different source distribution of primary electrons and primary positrons and nuclei
	std::vector<double> *parameters;
	int source_model;
	if (priElecStr != particle.name && priPosiStr != particle.name) {   //AWS20131023

		source_model = galdef.source_model;
		parameters  = &galdef.source_parameters;

	} else {

		if (priElecStr == particle.name) { //AWS20131023
			source_model = galdef.     source_model_electron;
			parameters  = &galdef.source_parameters_electron;
		}

		if (priPosiStr == particle.name) { //AWS20131023
			source_model = galdef.     source_model_positron;
			parameters  = &galdef.source_parameters_positron;
		}

	}


	if(galdef.n_spatial_dimensions==2) {


		/* ECC Feb 11 2015
		   Calculate the fraction of sources wihtin the CMZ
		*/

		double total_sources=0;
		double src_100=0;
		double src_300=0;
		double src_500=0;

		for(int ir=0; ir<particle.n_rgrid; ir++) {
			for(int iz=0; iz<particle.n_zgrid; iz++) {
				double sourceval = source_distribution(particle. r[ir], 0, particle.z[iz], source_model, *parameters);
				total_sources += sourceval*particle.r[ir];
				if (particle.r[ir] <= .1)
					src_100 += sourceval*particle.r[ir];
				if (particle.r[ir] <= .3)
					src_300 += sourceval*particle.r[ir];
				if (particle.r[ir] <= .5)
					src_500 += sourceval*particle.r[ir];
			}
		}

		cout << "============================================================" <<  endl;
		cout << "============================================================" <<  endl;
		cout << "Fraction of total sources within 100,300,500 pc of GC: "
		     << src_100/total_sources << ", "
		     << src_300/total_sources << ", "
		     << src_500/total_sources << endl;
		cout << "============================================================" <<  endl;
		cout << "============================================================" <<  endl;


		for(int ip=0; ip<particle.n_pgrid; ip++) {
			if(strcmp(galdef.inj_spectrum_type,"Etot")==0) spec_shape=pow(particle.Etot[ip],-g_2); // IMOS20000613
			else {                                                                                 // IMOS20000615
				/* AWS20131111
				 if(particle.rigidity[ip]< rigid_br0)                                   // IMOS20031012
				        spec_shape =pow(particle.rigidity[ip]/rigid_br0,-g_0) *pow(rigid_br0/rigid_br,-g_1);
				     if(rigid_br0<= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
				        spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_1);
				     if(rigid_br <= particle.rigidity[ip])
				        spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_2);
				*/

				if (particle.rigidity[ip] <  rigid_br0)
					spec_shape = pow(particle.rigidity[ip]/rigid_br0,-g_0)   * pow(rigid_br0/rigid_br, -g_1);

				if (particle.rigidity[ip] >= rigid_br0 && particle.rigidity[ip]< rigid_br )                     //AWS20131111
					spec_shape = pow(particle.rigidity[ip]/rigid_br,  -g_1);

				if (particle.rigidity[ip] >= rigid_br  && particle.rigidity[ip]< rigid_br2)
					spec_shape = pow(particle.rigidity[ip]/rigid_br,  -g_2);

				if (particle.rigidity[ip] >= rigid_br2)                                                         //AWS20131111
					spec_shape = pow(particle.rigidity[ip]/rigid_br2, -g_3)  * pow(rigid_br2/rigid_br, -g_2);




			}
//         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape*=particle.beta[ip];      // IMOS20000615
			if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape/=sqrt(1.+pow(2.e3/particle.rigidity[ip],2)); // IMOS20011210

			//Options to facilitate calculating the spectra piecewise for
			//gamma-ray and CR fitting
			if (particle.rigidity[ip] <= galdef.rigid_min)
				spec_shape = 0;
			if (particle.rigidity[ip] > galdef.rigid_max)
				spec_shape = 0;

			int ir=0, iz=particle.n_zgrid/2;
			if(galdef.source_specification==1) {
				particle.primary_source_function.d2[ir][iz].s[ip]=particle.primary_abundance*spec_shape;
				continue;
			}
			for(ir=0; ir<particle.n_rgrid; ir++) {
				if(galdef.source_specification==2) {
					particle.primary_source_function.d2[ir][iz].s[ip]=particle.primary_abundance*spec_shape;
					continue;
				}
				for(iz=0; iz<particle.n_zgrid; iz++) {
					if (0 == galdef.source_specification) {

						particle.primary_source_function.d2[ir][iz].s[ip] =
						    source_distribution(particle. r[ir], 0, particle.z[iz], source_model, *parameters)*particle.primary_abundance*spec_shape;



						//cout << particle.name << " " << ir << " " << iz << " " << ip << " " << particle.primary_source_function.d2[ir][iz].s[ip] << endl;

						if (galdef.verbose >= 1) {
							ost.str("");
							ost<<"r z source_distribution  "<<particle.r[ir]<<" "<<particle.z[iz]
							   <<" "<<source_distribution(     particle.r[ir],   0.0,particle.z[iz], source_model, *parameters);
							INFO(ost.str());
						}
					}
				}
			}
		}
	}

	//if (priElecStr == particle.name)
	//exit(0);

	if(galdef.n_spatial_dimensions==3) {

		for(int ip=0; ip<particle.n_pgrid; ip++) {
			if(strcmp(galdef.inj_spectrum_type,"Etot")==0) spec_shape=pow(particle.Etot[ip],-g_2); // IMOS20000613
			else {                                                                                 // IMOS20000615
				/* AWS20131111
				 if(particle.rigidity[ip]< rigid_br0)                                   // IMOS20031012
				        spec_shape =pow(particle.rigidity[ip]/rigid_br0,-g_0) *pow(rigid_br0/rigid_br,-g_1);

				     if(rigid_br0<= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
				        spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_1);

				     if(rigid_br <= particle.rigidity[ip])
				        spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_2);
				*/


				if (particle.rigidity[ip] <  rigid_br0)
					spec_shape = pow(particle.rigidity[ip]/rigid_br0,-g_0)   * pow(rigid_br0/rigid_br, -g_1);

				if (particle.rigidity[ip] >= rigid_br0 && particle.rigidity[ip]< rigid_br )                     //AWS20131111
					spec_shape = pow(particle.rigidity[ip]/rigid_br,  -g_1);

				if (particle.rigidity[ip] >= rigid_br  && particle.rigidity[ip]< rigid_br2)
					spec_shape = pow(particle.rigidity[ip]/rigid_br,  -g_2);

				if (particle.rigidity[ip] >= rigid_br2)                                                         //AWS20131111
					spec_shape = pow(particle.rigidity[ip]/rigid_br2, -g_3)  * pow(rigid_br2/rigid_br, -g_2);

			}
//         if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape*=particle.beta[ip];      // IMOS20000615
			if(strcmp(galdef.inj_spectrum_type,"beta_rig")==0) spec_shape/=sqrt(1.+pow(2.e3/particle.rigidity[ip],2)); // IMOS20011

			//Options to facilitate calculating the spectra piecewise for
			//gamma-ray and CR fitting
			if (particle.rigidity[ip] <= galdef.rigid_min)
				spec_shape = 0;
			if (particle.rigidity[ip] > galdef.rigid_max)
				spec_shape = 0;

			int ix=particle.n_xgrid/2, iy=particle.n_ygrid/2, iz=particle.n_zgrid/2;

			if(galdef.use_symmetry==1) ix=iy=iz=0;

			if(galdef.source_specification==1)
				particle.primary_source_function.d3[ix][iy][iz].s[ip]=particle.primary_abundance;

			double total_sources = 0;
			// ECC20141112: If the source distribution is multiple components, we need the relative normalization to be correct.
			// Integrate over the CO distribution and the parametrized version to count total number of sources in each.
			// Find total number of sources

			// Get grid spacing
			double dx = particle.x[1]-particle.x[0];
			double dy = particle.y[1]-particle.y[0];
			double dz = particle.z[1]-particle.z[0];

			if (galdef.spiral_fraction>0) {
				// if the normalizations have not been computed before, do it here.
				if ((source_norm_CO == -1) || (source_norm_param == -1) ) {
					source_norm_CO=0;
					source_norm_param=0;

					INFO("computing normalizations for non-zero spiral source fraction.");
					for(ix=0; ix<particle.n_xgrid; ix++) {
						for(iy=0; iy<particle.n_ygrid; iy++) {
							for(iz=0; iz<particle.n_zgrid; iz++) {

								// Average the gas density over the grid cell in a 10x10x10 array
								for(int ixx=0; ixx<3; ixx++) {
									for(int iyy=0; iyy<3; iyy++) {
										for(int izz=0; izz<3; izz++) {
											double x = particle.x[ix]+(ixx-1.)/3.*dx;
											double y = particle.y[iy]+(iyy-1.)/3.*dy;
											double z = particle.z[iz]+(izz-1.)/3.*dz;
											double r = pow(x*x+y*y, .5);

											double rho_H2 = interp_CO_gas_cube(x,y,z,2e20,COCubeData);
											if (rho_H2 > galdef.kennicutt_threshold) {
												source_norm_CO    += pow(rho_H2, galdef.kennicutt_index )/27;
											}

											source_norm_param += source_distribution(x,y,z,source_model,*parameters)/27;
										}//izz
									}//iyy
								}//ixx
							}//iz
						}//iy
					}//ix

					cout << "source_norm_CO, source_norm_param:" << source_norm_CO << " " << source_norm_param << endl;
					//cout << "kennicutt_index:" << galdef.kennicutt_index << " " << pow(interp_CO_gas_cube(0,0,0,1,COCubeData), galdef.kennicutt_index ) << endl;
				}// if
			}

			for(ix=0; ix<particle.n_xgrid; ix++) {
				for(iy=0; iy<particle.n_ygrid; iy++) {
					if(galdef.source_specification==2)
						particle.primary_source_function.d3[ix][iy][iz].s[ip]=particle.primary_abundance;

					for(iz=0; iz<particle.n_zgrid; iz++) {
						if(galdef.source_specification==0 && galdef.spiral_fraction==0)
							particle.primary_source_function.d3[ix][iy][iz].s[ip] =
							    source_distribution(particle.x[ix],particle.y[iy],particle.z[iz],
							                        source_model,*parameters)*particle.primary_abundance;
						// ECC20141112
						else if(galdef.source_specification==0 && galdef.spiral_fraction>0) {
							//if (spiral_fraction_particle_init == 0)
							//{
							for(int ixx=0; ixx<3; ixx++) {
								for(int iyy=0; iyy<3; iyy++) {
									for(int izz=0; izz<3; izz++) {
										double x = particle.x[ix]+(ixx-1.)/3.*dx;
										double y = particle.y[iy]+(iyy-1.)/3.*dy;
										double z = particle.z[iz]+(izz-1.)/3.*dz;
										double r = pow(x*x+y*y, .5);
										double rho_H2 = interp_CO_gas_cube(x,y,z,2e20,COCubeData);

										particle.primary_source_function.d3[ix][iy][iz].s[ip] +=
										    (1-galdef.spiral_fraction)/source_norm_param*source_distribution(x,y,z,source_model,*parameters)/27*particle.primary_abundance;

										if (rho_H2 > galdef.kennicutt_threshold) {

											// If we are in the CMZ, multiply by the CMZ_multiplier
											if ( (fabs(x)<.3) && (fabs(y)<.3) && (fabs(z)<.3)) {
												particle.primary_source_function.d3[ix][iy][iz].s[ip] +=
												    galdef.CMZ_multiplier*galdef.spiral_fraction/source_norm_CO*pow(rho_H2, galdef.kennicutt_index )/27*particle.primary_abundance;
											} else {
												particle.primary_source_function.d3[ix][iy][iz].s[ip] +=
												    galdef.spiral_fraction/source_norm_CO*pow(rho_H2, galdef.kennicutt_index )/27*particle.primary_abundance;
											}
										}
										//galdef.spiral_fraction/source_norm_CO*pow(interp_CO_gas_cube(x,y,z,fX_CO(r),COCubeData), galdef.kennicutt_index )/27*particle.primary_abundance;

										// spiral_fraction_particle.primary_source_function.d3[ix][iy][iz].s[ip] +=
										// (1-galdef.spiral_fraction)/source_norm_param*source_distribution(x,y,z,source_model,*parameters)/27*particle.primary_abundance;

										// spiral_fraction_particle.primary_source_function.d3[ix][iy][iz].s[ip] +=
										// galdef.spiral_fraction/source_norm_CO*pow(interp_CO_gas_cube(x,y,z,2e20,COCubeData), galdef.kennicutt_index )/27*particle.primary_abundance;
										//galdef.spiral_fraction/source_norm_CO*pow(interp_CO_gas_cube(x,y,z,fX_CO(r),COCubeData), galdef.kennicutt_index )/27*particle.primary_abundance;

									}//izz
								}//iyy
							}//ixx

							//}//end spiral_fraction_particle_init
							// else if the first particle has already been initialized, assume that the rest of them follow the same distributions
							// and multiply by the relative abundance w.r.t. gcr[0]
							// gcr[0] is the particle that was already initialized.  Just multiply the new species by the relative abundances
							//else{
							//   particle.primary_source_function.d3[ix][iy][iz].s[ip] = spiral_fraction_particle.primary_source_function.d3[ix][iy][iz].s[ip]*particle.primary_abundance/spiral_fraction_particle.primary_abundance;
							//}

						}


						// source spectral index dispersion

						double spec_dg_ratio;                             //AWS20010411
						spec_dg_ratio=1;                                  //AWS20080307

						if(galdef.SNR_events==1) {                         //AWS20080307
							if(strcmp(particle.name,"primary_electrons")==0) { //AWS20010411
								spec_dg_ratio=
								    pow(particle.rigidity[ip]/galdef.SNR_electron_dgpivot,1.*galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]);

								if(galdef.verbose==-501) { // selectable debug
									ost.str("");
									ost<<"SNR_electron_dg="<<galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]
									   <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio;
									INFO(ost.str());
								}
							}//electrons

							if(strcmp(particle.name,"primary_electrons")!=0) { //AWS20010411
								spec_dg_ratio=
								    pow(particle.rigidity[ip]/galdef.SNR_nuc_dgpivot,     1.*galaxy.SNR_nuc_dg.     d3[ix][iy][iz].s[0]);

								if(galdef.verbose==-501) { // selectable debug
									ost.str("");
									ost<<"SNR_nuc_dg="<<galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0]
									   <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio;
									INFO(ost.str());
								}
							}//nucleons
						}// if galdef.SNR_events

						// ECC 1/28/2016
						// If we are within 300 pc of the Galactic center, then use the inner injection spectrum specification
						// Otherwise use the usual injection spectrum.
						if ((particle.x[ix]*particle.x[ix] + particle.y[iy]*particle.y[iy]) < (.5*.5)
						    && ((priElecStr == particle.name) || strcmp(particle.name,"Hydrogen_1")==0 ) && (1==0)) { // Only set primary H and electron spec differently
							double spec_shape_inner;
							// if primary electrons
							if (priElecStr == particle.name) {
								if (particle.rigidity[ip] <  rigid_br0)
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br0,-g_0_inner) * pow(rigid_br0/rigid_br, -g_1_inner);

								if (particle.rigidity[ip] >= rigid_br0 && particle.rigidity[ip]< rigid_br )                     //AWS20131111
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br,  -g_1_inner);
								if (particle.rigidity[ip] >= rigid_br  && particle.rigidity[ip]< rigid_br2)
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br,  -g_2);
								if (particle.rigidity[ip] >= rigid_br2)                                                         //AWS20131111
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br2, -g_3)  * pow(rigid_br2/rigid_br, -g_2);
							} else {
								if (particle.rigidity[ip] <  rigid_br0)
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br0,-g_0_inner)   * pow(rigid_br0/rigid_br, -g_1_inner);
								if (particle.rigidity[ip] >= rigid_br0 && particle.rigidity[ip]< rigid_br )                     //AWS20131111
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br,  -g_1_inner);

								if (particle.rigidity[ip] >= rigid_br  && particle.rigidity[ip]< rigid_br2)
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br,  -g_2_inner);
								if (particle.rigidity[ip] >= rigid_br2)                                                         //AWS20131111
									spec_shape_inner = pow(particle.rigidity[ip]/rigid_br2, -g_3)  * pow(rigid_br2/rigid_br, -g_2_inner);
							}// end else

							//gamma-ray and CR fitting
							if (particle.rigidity[ip] <= galdef.rigid_min)
								spec_shape_inner = 0;
							if (particle.rigidity[ip] > galdef.rigid_max)
								spec_shape_inner = 0;

							particle.primary_source_function.d3[ix][iy][iz].s[ip]*=spec_shape_inner*spec_dg_ratio;
							// cout << particle.name << " Assigning Inner Source Spectrum.  rig/spec_shape_inner/spec_shape: " << particle.rigidity[ip] << " / " << spec_shape_inner << "/" << spec_shape << endl;
						}// end if electron or nucleus
						// if not in the inner region, keep things the same.
						else {
							particle.primary_source_function.d3[ix][iy][iz].s[ip]*=spec_shape*spec_dg_ratio;
						}// end else
					} //iz
				} //iy
			} //ix
		} //ip

	} // 3D

// CASE: PRIMARY NUCLEI                                                         AWS20000601.1

	if(strcmp(particle.species,"nucleus")==0) {                               // IMOS20000601.1
		particle.primary_source_function*= pow(particle.A, g_2-1);             // IMOS20000613.10
		if(strcmp(galdef.inj_spectrum_type,"Etot")!=0)                         // IMOS20000613.9
			particle.primary_source_function*= pow(fabs(1.*particle.Z),-g_2);   // AWS20000601.2
	}

	if(strcmp(particle.name,"primary_electrons")!=0 && strcmp(particle.name,"primary_positrons")!=0) { //AWS20131023
		particle.primary_source_function *= galdef.source_normalization;
		ost.str("");
		ost<<" >>>>>>>>>>>>>>>>>>> norm "<<galdef.source_normalization<<" >>>>>>>>>>>>>>>>>>>";
		INFO(ost.str());
	}

// CASE: PRIMARY ELECTRONS                                                      IMOS20031016

	if(strcmp(particle.name,"primary_electrons")==0) {
		particle.primary_source_function *= galdef.electron_source_normalization;
		ost.str("");
		ost<<" >>>>>>>>>> electron_norm "<<galdef.electron_source_normalization<<" >>>>>>>>>>>>>>>>>>>";
		INFO(ost.str());
	}

// CASE: PRIMARY POSITRONS                                                      AWS20131023

	if(strcmp(particle.name,"primary_positrons")==0) {
		particle.primary_source_function *= galdef.positron_source_normalization;
		ost.str("");
		ost<<" >>>>>>>>>> positron_norm "<<galdef.positron_source_normalization<<" >>>>>>>>>>>>>>>>>>>";
		INFO(ost.str());
	}


// ASSIGNING DIFFUSION COEFFICIENT

	INFO("      assigning diffusion coefficient");

// compute beta at break rigidity so that formula will give galdef.D0_xx at this point
	double Ekin_br=-1.;                             // so that p will be used in kinematic
	double p_br=galdef.D_rigid_br*fabs(1.*particle.Z);
	double Etot_br, beta_br, gamma_br, rigidity_br; // output of kinematic
	char species[10];
////////////////////////V//IMOS20030214 all region
	int iprotons=-1;

	if(galdef.diff_reacc > 5) {
// identify CR protons
		for(int i=0; i<n_species; i++)
			if(101==100*gcr[i].Z+gcr[i].A) {
				iprotons=i;
				ost.str("");
				ost<<"  CR protons found as species #"<<iprotons;
				INFO(ost.str());
				break;
			}
		if(iprotons==-1) {
			WARNING("CR protons not found!");
			return 1;
		}
		ost.str("");
		ost<<"create_transport_arrays>> "<<particle.Z*100+particle.A<<" "<<particle.p[0];
		INFO(ost.str());

// Zero approximation proton spectrum for calculation of damping
		if(gcr[iprotons].cr_density.max() == 0) {
			if(galdef.n_spatial_dimensions==2)
				for(int ir=0; ir<particle.n_rgrid; ir++)
					for(int iz=0; iz<particle.n_zgrid; iz++)
						for(int ip=0; ip<particle.n_pgrid; ip++) {
							gcr[iprotons].cr_density.d2[ir]    [iz].s[ip] =
							    (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
							    *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
//cout<<"create_transport_arrays>> "<<gcr[iprotons].cr_density.d2[ir]    [iz].s[ip]<<endl;
						}

			if(galdef.n_spatial_dimensions==3)
				for(int ix=0; ix<particle.n_xgrid; ix++)
					for(int iy=0; iy<particle.n_ygrid; iy++)
						for(int iz=0; iz<particle.n_zgrid; iz++)
							for(int ip=0; ip<particle.n_pgrid; ip++)
								gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip] =
								    (1.-(galdef.z_min+iz*galdef.dz)/galdef.z_max)
								    *galdef.proton_norm_flux *pow(gcr[iprotons].Etot[ip]/galdef.proton_norm_Ekin, -2.75);
		}
		key=1;
	}
//   cout<<"create_transport_arrays>> "<<iprotons<<" "<<gcr[iprotons].cr_density.d2[9][9].s[10]<< " "<<protons.d2[9][9].s[10]<<endl;
//   for(int ip=0; ip<particle.n_pgrid; ip++) cout<<" "<<protons.d2[9][9].s[ip];
//   cout<<endl;
////////////////////////^//IMOS20030214

	strcpy(species,"nucleus");
	if(particle.A==0) strcpy(species,"electron");

// IMOS20000810.2
	if(kinematic(particle.Z,particle.A,species,p_br,Ekin_br,Etot_br,beta_br,gamma_br,rigidity_br,0)) exit(1);

	ost.str("");
	ost<<" beta at break rigidity="<<beta_br<<"   rigidity_br="<<rigidity_br
	   <<"?= galdef.D_rigid_br="<<galdef.D_rigid_br;
	INFO(ost.str());

	beta_br=1.0;  // to simulate fortran implementation

	if(galdef.n_spatial_dimensions==2) {
		ost.str("");
		ost<<" Calculating diffusion coefficients (ir, iz, ip): ("<<particle.n_rgrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
		INFO(ost.str());
		#pragma omp parallel for schedule(dynamic) default(shared)
		for(int ir=0; ir<particle.n_rgrid; ir++)
			for(int iz=0; iz<particle.n_zgrid; iz++)
				for(int ip=particle.n_pgrid-1; ip>=0; ip--) { // IMOS20060330 changed to reverse order (Wave-particle interactions)
					D_xx(particle,iprotons,ir, 0, 0,iz,ip); //array assigned in D_xx IMOS20030129
					particle.Dpp.d2[ir][iz].s[ip] =D_pp(particle.p[ip],galdef.D_g_1,galdef.v_Alfven,particle.Dxx.d2[ir][iz].s[ip]);
				}
	}
	if(galdef.n_spatial_dimensions==3) {
		ost.str("");
		ost<<" Calculating diffusion coefficients (ix, iy, iz, ip): ("<<particle.n_xgrid<<", "<<particle.n_ygrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
		INFO(ost.str());
		#pragma omp parallel for schedule(dynamic) default(shared)
		for(int ix=0; ix<particle.n_xgrid; ix++)
			for(int iy=0; iy<particle.n_ygrid; iy++)
				for(int iz=0; iz<particle.n_zgrid; iz++)
					for(int ip=particle.n_pgrid-1; ip>=0; ip--) { // IMOS20060330 changed to reverse order
						D_xx(particle,iprotons, 0,ix,iy,iz,ip); //array assigned in D_xx IMOS20030129
						particle.Dpp.d3[ix][iy][iz].s[ip] =D_pp(particle.p[ip], galdef.D_g_1, galdef.v_Alfven, particle.Dxx.d3[ix][iy][iz].s[ip]);
					}  // p
	}
	if(galdef.verbose>=2) {
		ost.str("");
		ost<< "spatial   diffusion coefficient Dxx  for species "<<particle.name;
		INFO(ost.str());
		particle.Dxx.print();
		ost.str("");
		ost<< "momentum diffusion coefficient Dpp  for species "<<particle.name;
		INFO(ost.str());
		particle.Dpp.print();
	}

	// ASSIGNING CONVECTION ARRAYS                                         AWS20131008

	if(galdef.n_spatial_dimensions==2) {
		ost.str("");
		ost<<" Calculating convection velocity    (ir, iz, ip): ("<<particle.n_rgrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
		INFO(ost.str());

		for(int ir=0; ir<particle.n_rgrid; ++ir)
			for(int iz=0; iz<particle.n_zgrid; ++iz)
				for(int ip=particle.n_pgrid-1; ip>=0; --ip) {
					if(galdef.convection==1) { //AWS20130320
						particle.   v_conv.d2[ir][iz].s[ip] = galdef.v0_conv + galdef.dvdz_conv  *   particle.z[iz] ; // original definition
						particle.dvdz_conv.d2[ir][iz].s[ip] =                  galdef.dvdz_conv;                      //AWS20130502
					}

					if(galdef.convection==2) { // linear for |z|>z0, zero for |z|<z0  AWS20130320
						particle.   v_conv.d2[ir][iz].s[ip] =  0.;
						if(particle.z[iz] >=  galdef.z0_conv)      particle.   v_conv.d2[ir][iz].s[ip] =  galdef.v0_conv + galdef.dvdz_conv * (abs(particle.z[iz]) - galdef.z0_conv) ;
						if(particle.z[iz] <= -galdef.z0_conv)      particle.   v_conv.d2[ir][iz].s[ip] = -galdef.v0_conv - galdef.dvdz_conv * (abs(particle.z[iz]) - galdef.z0_conv) ;

						particle.dvdz_conv.d2[ir][iz].s[ip] =  0.;                                    //AWS20130502
						if(abs(particle.z[iz]) >=  galdef.z0_conv) particle.dvdz_conv.d2[ir][iz].s[ip] =  galdef.dvdz_conv;                      //AWS20130502

					}

					if(galdef.convection==3) { // smooth increase from zero at z=0 to v0_conv at large z. 0.5*v0_conv at z0_conv  AWS20130429
						particle.v_conv.d2[ir][iz].s[ip] =  galdef.v0_conv * 0.5 * (1.0 + tanh( 4.0*abs(particle.z[iz]/galdef.z0_conv) - 4.0) );
						if(particle.z[iz] < 0.)   particle.v_conv.d2[ir][iz].s[ip]*= -1.; // velocity negative below plane

						particle.dvdz_conv.d2[ir][iz].s[ip] =
						    galdef.v0_conv *0.5 *4.0  /galdef.z0_conv *pow(cosh(4.0*abs(particle.z[iz]/galdef.z0_conv) - 4.0) ,-2);

					}

					if(galdef.convection==4) { // constant radial (3D) wind until fermi dirac style shutoff

						// 3D radius
						double r3 = sqrt(particle.r[ir]*particle.r[ir] + particle.z[iz]*particle.z[iz]);
						// Shut off wind with fermi dirac distribution with stall transition-width 200 pc and and 2 kpc stall height
						double vr = 0;
						if (r3<4.) {
							vr = galdef.v0_conv * 1/(exp((r3-2.)/.2)+1);
						}
						if (isnan(vr)) {
							//cout << "x,y,z ---- r3, vr : " << particle.x[ix] << ", " << particle.y[iy] << ", " << particle.z[iz] << " -------- " << r3 << ", " << vr << endl;
						}

						if (r3!=0) {
							// Assign v_conv_i
							particle.v_conv_x.d2[ir][iz].s[ip] =  particle.r[ir]/r3*vr; // x is used for the r direction
							particle.v_conv_z.d2[ir][iz].s[ip] =  particle.z[iz]/r3*vr; //
						} else {
							// Assign v_conv_i
							// particle.v_conv_x.d3[ix][iy][iz].s[ip] =  vr;
							// particle.v_conv_y.d3[ix][iy][iz].s[ip] =  vr;
							// particle.v_conv_z.d3[ix][iy][iz].s[ip] =  vr;
							// What do we do exactly at the origin?  For now just set to zero..
							particle.v_conv_x.d2[ir][iz].s[ip] =  0;
							particle.v_conv_z.d2[ir][iz].s[ip] =  0;
						}
					}
				}
	}

	if(galdef.n_spatial_dimensions==3) {
		ost.str("");
		ost<<" Calculating convection velocity    (ix, iy, iz, ip): ("<<particle.n_xgrid<<", "<<particle.n_ygrid<<", "<<particle.n_zgrid<<", "<<particle.n_pgrid<<")";
		INFO(ost.str());

		for(int ix=0; ix<particle.n_xgrid; ++ix)
			for(int iy=0; iy<particle.n_ygrid; ++iy)
				for(int iz=0; iz<particle.n_zgrid; ++iz)
					for(int ip=particle.n_pgrid-1; ip>=0; --ip) {
						if(galdef.convection==1) { //AWS20130320
							particle.   v_conv.d3[ix][iy][iz].s[ip] = galdef.v0_conv + galdef.dvdz_conv  *   particle.z[iz] ; // original definition
							particle.dvdz_conv.d3[ix][iy][iz].s[ip] =                  galdef.dvdz_conv;                      //AWS20130502
						}

						if(galdef.convection==2) { // linear for |z|>z0, zero for |z|<z0  AWS20130320
							particle.v_conv.d3[ix][iy][iz].s[ip] =  0.;
							if(particle.z[iz] >=  galdef.z0_conv)  particle.v_conv.d3[ix][iy][iz].s[ip] =  galdef.v0_conv + galdef.dvdz_conv * (abs(particle.z[iz]) - galdef.z0_conv) ;
							if(particle.z[iz] <= -galdef.z0_conv)  particle.v_conv.d3[ix][iy][iz].s[ip] = -galdef.v0_conv - galdef.dvdz_conv * (abs(particle.z[iz]) - galdef.z0_conv) ;

							particle.dvdz_conv.d3[ix][iy][iz].s[ip] =  0.;                                    //AWS20130502
							if(abs(particle.z[iz]) >=  galdef.z0_conv) particle.dvdz_conv.d3[ix][iy][iz].s[ip] =  galdef.dvdz_conv;                      //AWS20130502
						}

						if(galdef.convection==3) { // smooth increase from zero at z=0 to v0_conv at large z. 0.5*v0_conv at z0_conv  AWS20130429
							particle.v_conv.d3[ix][iy][iz].s[ip] =  galdef.v0_conv * 0.5 * (1.0 + tanh( 4.0*abs(particle.z[iz]/galdef.z0_conv) - 4.0) );
							if(particle.z[iz] < 0.)   particle.v_conv.d3[ix][iy][iz].s[ip]*= -1.; // velocity negative below plane


							particle.dvdz_conv.d3[ix][iy][iz].s[ip] =
							    galdef.v0_conv *0.5 *4.0  /galdef.z0_conv *pow(cosh(4.0*abs(particle.z[iz]/galdef.z0_conv) - 4.0) ,-2);

						}
						// ECC 1.28.2016
						if(galdef.convection==4) { // constant radial (3D) wind until fermi dirac style shutoff
							// 3D radius
							double r3 = sqrt(particle.x[ix]*particle.x[ix] + particle.y[iy]*particle.y[iy] + particle.z[iz]*particle.z[iz]);
							// Shut off wind with fermi dirac distribution with stall transition-width 200 pc and and 2 kpc stall height
							double vr = 0;
							if (r3<4.) {
								vr = galdef.v0_conv * 1/(exp((r3-2.)/.2)+1);
							}
							if (isnan(vr)) {
								cout << "x,y,z ---- r3, vr : " << particle.x[ix] << ", " << particle.y[iy] << ", " << particle.z[iz] << " -------- " << r3 << ", " << vr << endl;
							}

							if (r3!=0) {
								// Assign v_conv_i
								particle.v_conv_x.d3[ix][iy][iz].s[ip] =  particle.x[ix]/r3*vr;
								particle.v_conv_y.d3[ix][iy][iz].s[ip] =  particle.y[iy]/r3*vr;
								particle.v_conv_z.d3[ix][iy][iz].s[ip] =  particle.z[iz]/r3*vr;
							} else {
								// Assign v_conv_i
								// particle.v_conv_x.d3[ix][iy][iz].s[ip] =  vr;
								// particle.v_conv_y.d3[ix][iy][iz].s[ip] =  vr;
								// particle.v_conv_z.d3[ix][iy][iz].s[ip] =  vr;
								// What do we do exactly at the origin?  For now just set to zero..
								particle.v_conv_x.d3[ix][iy][iz].s[ip] =  0;
								particle.v_conv_y.d3[ix][iy][iz].s[ip] =  0;
								particle.v_conv_z.d3[ix][iy][iz].s[ip] =  0;
							}
						}

					}  // p
	}


	if(galdef.verbose==-504) { // selectable debug
		ost.str("");
		ost<< " convection velocity  for species "<<particle.name;
		INFO(ost.str());
		cout<< " convection velocity  for species "<<particle.name<<endl;;
		particle.v_conv.print();

		ost.str("");
		ost<< " convection dvdz      for species "<<particle.name;
		INFO(ost.str());
		cout<< " convection dvdz      for species "<<particle.name<<endl;
		particle.dvdz_conv.print();


	}








// ASSIGNING FRAGMENTATION RATE

	INFO("======== assigning fragmentation rate ======== ");
	int ZH=1, ZHe=2;                                                    //  IMOS20010816
	double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;          //  IMOS20010816

	particle.fragment=0.0;
	if(particle.A!=0) {
		double CSratio,CStot_ratio;

		for(int ip=0; ip<particle.n_pgrid; ip++) {
			// IMOS20000607 whole segment
			A1 = 1;                                                    // nucleus
			Z2 = particle.Z;
			A2 = particle.A;                         // - " -
			if(101==100*Z2+A2) {
				Z2 = 2;     // protons
				A2 = 4;
			}
			if(-99==100*Z2+A2) {
				A1 =-1;     // antiprotons
				Z2 = 2;
				A2 = 4;
			}
			nucleon_cs(galdef.total_cross_section,particle.Ekin[ip]*1.e-3,A1,Z2,A2, // AWS20010620
			           &PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);
			He_to_H_CS(particle.Ekin[ip]*1.e-3,particle.Z,particle.A,999,999,&CSratio,&CStot_ratio);

			fragment_p = PA_inel;                                      // nuclei
			fragment_He= PA_inel*CStot_ratio;                          // -"-

// ELECTRON ATTACHMENT/STRIPPING CROSS SECTION                                  IMOS20010816

			if(galdef.K_capture) {
				for(K_electron=0; K_electron<=galdef.K_capture; K_electron++)
					nucdata(galdef_network_par,particle.Z,particle.A,K_electron,particle.Z,particle.A,&Z2,&A2,&t_half[K_electron]);

				if(t_half[0] != t_half[1]) {
					Kcapture_cs(particle.Ekin[ip],particle.Z,ZH, &attach_H ,&strip_H );
					Kcapture_cs(particle.Ekin[ip],particle.Z,ZHe,&attach_He,&strip_He);
					if(particle.K_electron) {
						fragment_p += strip_H ;
						fragment_He+= strip_He;
					} else {
						fragment_p += attach_H ;
						fragment_He+= attach_He;
					}

					if(galdef.verbose==-502) { // selectable debug
						ost.str("");
						ost<<"create_transport_arrays: Ekin Z,A,K_electron,attach_H,strip_H: "
						   <<particle.Ekin[ip]<<" "<<particle.Z<<" "<<particle.A<<" "
						   <<particle.K_electron<<" "<<strip_H<<" "<<attach_H;
						INFO(ost.str());
					}
				}
			}
			if(101==100*particle.Z+particle.A) {                       // protons
				fragment_p = PP_inel;
				fragment_He= PA_inel;
			}
			if(-99==100*particle.Z+particle.A) {                       // antiprotons
				fragment_p = aPP_non+aPP_ann;
				fragment_He= aPA_non+aPA_ann;
			}

			if(galdef.n_spatial_dimensions==2)
				for(int ir=0; ir<particle.n_rgrid; ir++)                   // IMOS20010816
					for(int iz=0; iz<particle.n_zgrid; iz++)
						particle.fragment.d2[ir] [iz].s[ip]= particle.beta[ip]*C
						                                     *(galaxy.n_HI.d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0]+galaxy.n_HII.d2[ir] [iz].s[0])
						                                     *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;

			if(galdef.n_spatial_dimensions==3)
				for(int ix=0; ix<particle.n_xgrid; ix++)                // IMOS20010816
					for(int iy=0; iy<particle.n_ygrid; iy++)
						for(int iz=0; iz<particle.n_zgrid; iz++)
							particle.fragment.d3[ix][iy][iz].s[ip]= particle.beta[ip]*C
							                                        *(galaxy.n_HI.d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
							                                        *(fragment_p+galdef.He_H_ratio*fragment_He) *1.0e-27;
		}  //  p
	}  //  A!=0

	if(galdef.verbose>=2) {
		ost.str("");
		ost<< "fragmentation for species "<<particle.name;
		INFO(ost.str());
		particle.fragment.print();
	}

// ASSIGNING MOMENTUM LOSS RATE

	INFO("======== assigning momentum loss rate ======== ");

	if(galdef.n_spatial_dimensions==2) {
		for(int ir=0; ir<particle.n_rgrid; ir++) {
			for(int iz=0; iz<particle.n_zgrid; iz++) {
				for(int ip=0; ip<particle.n_pgrid; ip++) {
					double aion,coul;                                      // NUCLEONS

					if(particle.A!=0) particle.dpdt.d2[ir] [iz].s[ip]=
						    nucleon_loss(particle.Z,particle.A,particle.Etot[ip],
						                 galaxy.n_HI .d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0],
						                 galaxy.n_HII.d2[ir] [iz].s[0], galdef.He_H_ratio,
						                 &aion,  &coul) / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1

					if(particle.A==0) {
						double uevcm3=0., bevcm3, brem1,brem2,sync,cmptn;   // ELECTRONS
						// IMOS200008016
// test of electron propagation vs analytical calculations (only IC losses) IMOS20061030
						if(abs(galdef.DM_int0)==99) {
							particle.dpdt.  d2[ir] [iz].s[ip]=
							    electron_loss( particle.Etot[ip], 0., 0., galdef.He_H_ratio, galdef.DM_double7, 0.,
							                   &aion, &coul,&brem1,&brem2,&sync,&cmptn) *1.0e-6; // energy eV s-1 -> momentum MeV s-1
							continue;
						}
// end of the test area

						bevcm3=pow(galaxy.B_field.d2[ir][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3
						particle.dpdt.  d2[ir] [iz].s[ip]= electron_loss( particle.Etot[ip],
						                                   galaxy.n_HI .d2[ir] [iz].s[0] +2*galaxy.n_H2.d2[ir] [iz].s[0],
						                                   galaxy.n_HII.d2[ir] [iz].s[0], galdef.He_H_ratio, uevcm3, bevcm3,
						                                   &aion, &coul,&brem1,&brem2,&sync,&cmptn)
						                                   / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1
					}  //  A==0
// cout<<" dpdt="<<particle.dpdt.d2[ix]    [iz].s[ip]<<" aion="<<aion<<endl;
				}  //  p
			}  //  z
		}  //  r
	}

	if(galdef.n_spatial_dimensions==3) {
		for(int ix=0; ix<particle.n_xgrid; ix++) {
			for(int iy=0; iy<particle.n_ygrid; iy++) {
				for(int iz=0; iz<particle.n_zgrid; iz++) {
					for(int ip=0; ip<particle.n_pgrid; ip++) {
						double aion,coul;                                        // NUCLEONS

						if(particle.A!=0) particle.dpdt.d3[ix][iy][iz].s[ip]=
							    nucleon_loss(particle.Z,particle.A,particle.Etot[ip],
							                 galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
							                 galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio,
							                 &aion,&coul) / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1

						if(particle.A==0) {
							double uevcm3=0., bevcm3=0.,brem1,brem2,sync,cmptn;   // ELECTRONS
							// IMOS200008016
							bevcm3=pow(galaxy.B_field.d3[ix][iy][iz].s[0]*10000.,2)/8./Pi *erg_to_eV;// mag. energy density eV cm-3
							particle.  dpdt.d3[ix][iy][iz].s[ip]= electron_loss(particle.Etot[ip],
							                                      galaxy.n_HI .d3[ix][iy][iz].s[0] +2*galaxy.n_H2.d3[ix][iy][iz].s[0],
							                                      galaxy.n_HII.d3[ix][iy][iz].s[0],galdef.He_H_ratio, uevcm3, bevcm3,
							                                      &aion,&coul,&brem1,&brem2,&sync,&cmptn)
							                                      / particle.beta[ip]*1.0e-6; // energy eV s-1 -> momentum MeV s-1
						}  //  A==0
//  cout<<" dpdt="<<particle.dpdt.d3[ix][iy][iz].s[ip]<<" p="<<particle.p[ip] <<endl;
					}  //  p
				}  //  z
			}  //  y
		}  //  x
	}

// IF ELECTRON ADD KLEIN_NISHINA LOSSES

	if(abs(galdef.DM_int0)!=99) if(particle.A==0) e_KN_loss(particle);  // MeV s-1 IMOS20061030

	if(galdef.verbose>=2) {
		ost.str("");
		ost<< "dpdt for species "<<particle.name;
		INFO(ost.str());
		particle.dpdt.print();
	}

// ASSIGNING DECAY RATE

	if(particle.t_half!=0.0) {
		INFO("======== assigning decay rate ======== ");

		if(galdef.n_spatial_dimensions==2) {
			for(int ir=0; ir<particle.n_rgrid; ir++) {
				for(int iz=0; iz<particle.n_zgrid; iz++) {
					for(int ip=0; ip<particle.n_pgrid; ip++)
						particle.decay.d2[ir][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
				}  // z
			}  //  r
		}

		if(galdef.n_spatial_dimensions==3) {
			for(int ix=0; ix<particle.n_xgrid; ix++) {
				for(int iy=0; iy<particle.n_ygrid; iy++) {
					for(int iz=0; iz<particle.n_zgrid; iz++) {
						for(int ip=0; ip<particle.n_pgrid; ip++)
							particle.decay.d3[ix][iy][iz].s[ip]=1.0/(particle.gamma[ip]*particle.t_half*year2sec/log(2.0));
					}  //  z
				}  //  y
			}  //  x
		}

		if(galdef.verbose>=1) {
			ost.str("");
			ost<< "decay for species "<<particle.name;
			INFO(ost.str());
			particle.decay.print();
		}
	}  //  t_half!=0.0

	if(galdef.verbose>=1) {
		ost.str("");
		ost<< "primary source function for species "<<particle.name;
		INFO(ost.str());
		particle.primary_source_function.print();
	}

//particle.print();
	ost.str("");
	ost<<"============== completed creation of transport arrays for "<<particle.name;
	INFO(ost.str());
	INFO("Exit");
	return stat;
}




