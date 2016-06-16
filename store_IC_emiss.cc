
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_emiss.cc *                       galprop package * 4/14/2000 
// * Modified to output 3D BJB20160615
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"
#include "fitsio.h" 

#include <ErrorLogger.h>

void Galprop::store_IC_emiss(const string& type) {

  INFO("Entry");

  fitsfile* fptr = 0;

  const long naxis = 5;

  long naxes[naxis];

  assert(galaxy.n_spatial_dimensions == 2 || galaxy.n_spatial_dimensions == 3);

  if(galaxy.n_spatial_dimensions==2) {
  	naxes[0]=galaxy.n_rgrid;
  	naxes[1]=galaxy.n_zgrid;
    naxes[2]=galaxy.n_E_gammagrid;
  	naxes[3]=galaxy.n_ISRF_components + 1;
  	naxes[4]=1;
  }
  
  if(galaxy.n_spatial_dimensions==3) {
  	naxes[0]=galaxy.n_xgrid;
  	naxes[1]=galaxy.n_ygrid;             
  	naxes[2]=galaxy.n_zgrid;
  	naxes[3]=galaxy.n_E_gammagrid;
  	naxes[4]=galaxy.n_ISRF_components + 1;
  }
  
	/*axes[0] = (galaxy.n_spatial_dimensions == 2 ? galaxy.n_rgrid : galaxy.n_xgrid);
  axes[1] = galaxy.n_zgrid;
  axes[2] = galaxy.n_E_gammagrid;
  axes[3] = galaxy.n_ISRF_components + 1;*/

  const long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3]*naxes[4];

  valarray<float> array(0., nElements);

  const string filename = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "IC_emiss_" + galdef.galdef_ID + ".gz";

  int status = 0;

  fits_create_file(&fptr, filename.c_str(), &status);

  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

  // Write some keywords giving information about whether its an isotropic 
  // emissivity, anisotropic calculated for particular viewing location, etc.

  long isotropicIC = (type == "isotropic" ? 1 : 0);
  
  fits_update_key(fptr, TLONG, "ISOTROPIC", &isotropicIC, "Isotropic/Anisotropic cross section", &status);

  // for 3D case store x-dimension at y = 0 -- for now consistent with brem
  // emissivity storage etc. In future will export full 3D 

  int i=0;

	if(galaxy.n_spatial_dimensions==2) {
		for (int icomp = 0; icomp < naxes[3]; ++icomp)
  	for (int ip=0; ip<naxes[2]; ip++)
  	for (int iz=0; iz<naxes[1]; iz++)
		for (int ir=0; ir<naxes[0]; ir++) 
		{
			if (icomp < naxes[3]-1) {
				if (type=="isotropic")
					array[i]=galaxy.IC_iso_emiss[icomp].d2[ir][iz].s[ip];
				else
					array[i]=galaxy.IC_aniso_emiss->d2[ir][iz].s[ip];
				i++;
			}
			else {
				double sum=0;
				for (int ic=0; ic < naxes[3]-1; ++ic)
					if (type=="isotropic")
						sum+=galaxy.IC_iso_emiss[ic].d2[ir][iz].s[ip];
					else
						sum+=galaxy.IC_aniso_emiss->d2[ir][iz].s[ip];
				array[i]=sum;
			}
			array[i]*=galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
			i++;
		}
	}
 
	if(galaxy.n_spatial_dimensions==3) {
		for (int icomp = 0; icomp < naxes[4]; ++icomp)
		for (int ip=0; ip<naxes[3]; ip++)
		for (int iz=0; iz<naxes[2]; iz++)
		for (int iy=0; iy<naxes[1]; iy++)
		for (int ix=0; ix<naxes[0]; ix++)
		{
			if (icomp < naxes[3]-1) {
				if (type=="isotropic")
					array[i]=galaxy.IC_iso_emiss[icomp].d3[ix][iy][iz].s[ip];
				else
					array[i]=galaxy.IC_aniso_emiss->d3[ix][iy][iz].s[ip];
				i++;
			}
			else {
				double sum=0;
				for (int ic=0; ic < naxes[3]-1; ++ic)
					if (type=="isotropic")
						sum+=galaxy.IC_iso_emiss[ic].d3[ix][iy][iz].s[ip];
					else
						sum+=galaxy.IC_aniso_emiss->d3[ix][iy][iz].s[ip];
				array[i]=sum;
			}
			array[i]*=galaxy.E_gamma[ip]*galaxy.E_gamma[ip];
			i++;
		}
	}
	

	/*if(galaxy.n_spatial_dimensions==2) {
		for (unsigned int iComp = 0; iComp < naxes[3]; ++iComp) 
		for (unsigned int iP = 0; iP < naxes[2]; ++iP)
		for (unsigned int iZ = 0; iZ < naxes[1]; ++iZ)
		for (unsigned int iR = 0; iR < naxes[0]; ++iR) {

			unsigned int index = iComp*naxes[2]*naxes[1]*naxes[0] + iP*naxes[1]*naxes[0] + iZ*naxes[0] + iR; // explicit encoding -- use it!!

			if (iComp < axes[3] - 1) { // Individual components

			  if (type == "isotropic")
			    array[index] = (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_iso_emiss[iComp].d2[iR][iZ].s[iP] : galaxy.IC_iso_emiss[iComp].d3[iR][0][iZ].s[iP]);
			  else
			    array[index] = (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP] : galaxy.IC_aniso_emiss->d3[iR][0][iZ].s[iP]);

			} 
			else { // Total emissivity

				double sum = 0;

			  for (unsigned int iC = 0; iC < axes[3]-1; ++iC)
			    if (type == "isotropic")
						sum += (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_iso_emiss[iC].d2[iR][iZ].s[iP] : galaxy.IC_iso_emiss[iC].d3[iR][0][iZ].s[iP]);
			    else
						sum += (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP] : galaxy.IC_aniso_emiss->d3[iR][0][iZ].s[iP]);

			  array[index] = sum;

			}

			array[index] *= galaxy.E_gamma[iP]*galaxy.E_gamma[iP];

		}
	} */
  
  // Write the array of floats to the image 
  long fPixel = 1;
  fits_write_img(fptr, TFLOAT, fPixel, nElements, &array[0], &status);
  
  // write basic FITS keywords

  double crval1,crval2,crval3,crval4,crval5;
  double cdelt1,cdelt2,cdelt3,cdelt4,cdelt5;

  if(galaxy.n_spatial_dimensions==2) {
		crval1=galaxy.r_min;
		crval2=galaxy.z_min;
		crval3=log10(galaxy.E_gamma_min);
		crval4=1;
		crval5=1;
		
		cdelt1=galaxy.dr;
		cdelt2=galaxy.dz;
		cdelt3=log10(galaxy.E_gamma_factor);
		cdelt4=1;
		cdelt5=1;
	}
		
  if(galaxy.n_spatial_dimensions==3) {
		crval1=galaxy.x_min;
		crval2=galaxy.y_min;
		crval3=galaxy.z_min;
		crval4=log10(galaxy.E_gamma_min);
		crval5=1;
  
		cdelt1=galaxy.dx;
		cdelt2=galaxy.dy;
		cdelt3=galaxy.dz;
		cdelt4=log10(galaxy.E_gamma_factor);
		cdelt5=1;
	}
  
  /*double crval1 = (galaxy.n_spatial_dimensions == 2 ? galaxy.r_min : galaxy.x_min);
  double crval2 = galaxy.z_min;
  double crval3 = log10(galaxy.E_gamma_min);
  double crval4 = 1;

  double cdelt1 = (galaxy.n_spatial_dimensions == 2 ? galaxy.dr : galaxy.dx);
  double cdelt2 = galaxy.dz;
  double cdelt3 = log10(galaxy.E_gamma_factor);
  double cdelt4 = 1; */

  if(galaxy.n_spatial_dimensions==2) {
		fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of radial dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of Z dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of log10(energy grid/MeV)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of component numbering", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL5", &crval5,"Nothing", &status);
		
		fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of radial dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of Z dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of log10(energy grid/MeV)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of components", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT5", &cdelt5,"Nothing", &status);
  }
  
  if(galaxy.n_spatial_dimensions==3) {
		fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of X dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of Y dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of Z dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of log10(energy grid/MeV)", &status);
		fits_update_key(fptr, TDOUBLE, "CRVAL5", &crval5,"Start of component numbering", &status);
		
		fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of X dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of Y dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of Z dimension (kpc)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of log10(energy grid/MeV)", &status);
		fits_update_key(fptr, TDOUBLE, "CDELT5", &cdelt5,"Increment of components", &status);
  }
  
  fits_close_file(fptr, &status);   
  fits_report_error(stderr, status);

  INFO("Exit");

}
