
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_emiss.cc *                       galprop package * 4/14/2000
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

void Galprop::store_IC_emiss(const string& type)
{

	INFO("Entry");

	fitsfile* fptr = 0;

	assert(galaxy.n_spatial_dimensions == 2 || galaxy.n_spatial_dimensions == 3);

	const long nAxes = 5;
	long axes[nAxes];

	if(2==galaxy.n_spatial_dimensions) {
		//const long nAxes = 4;
		//long axes[nAxes];

		axes[0] = galaxy.n_rgrid;
		axes[1] = galaxy.n_zgrid;
		axes[2] = galaxy.n_E_gammagrid;
		axes[3] = galaxy.n_ISRF_components + 1;
		axes[4] = 1;

		//const long nElements = axes[0]*axes[1]*axes[2]*axes[3];
	} //2D

	if(3==galaxy.n_spatial_dimensions) {
		//const long nAxes = 5;
		//long axes[nAxes];

		axes[0] = galaxy.n_xgrid;
		axes[1] = galaxy.n_ygrid;
		axes[2] = galaxy.n_zgrid;
		axes[3] = galaxy.n_E_gammagrid;
		axes[4] = galaxy.n_ISRF_components + 1;

		//const long nElements = axes[0]*axes[1]*axes[2]*axes[3]*axes[4];
	} //3D

	const long nElements = axes[0]*axes[1]*axes[2]*axes[3]*axes[4];

	valarray<float> array(0., nElements);

	const string filename = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "IC_emiss_" + galdef.galdef_ID + ".gz";

	int status = 0;

	fits_create_file(&fptr, filename.c_str(), &status);

	fits_create_img(fptr, FLOAT_IMG, nAxes, axes, &status);

	// Write some keywords giving information about whether its an isotropic
	// emissivity, anisotropic calculated for particular viewing location, etc.

	long isotropicIC = (type == "isotropic" ? 1 : 0);

	fits_update_key(fptr, TLONG, "ISOTROPIC", &isotropicIC, "Isotropic/Anisotropic cross section", &status);

	if(2==galaxy.n_spatial_dimensions) {
		for (unsigned int iComp = 0; iComp < axes[3]; ++iComp) {
			for (unsigned int iP = 0; iP < axes[2]; ++iP) {
				for (unsigned int iZ = 0; iZ < axes[1]; ++iZ) {
					for (unsigned int iR = 0; iR < axes[0]; ++iR) {
						unsigned int index = iComp*axes[2]*axes[1]*axes[0] + iP*axes[1]*axes[0] + iZ*axes[0] + iR; // explicit encoding -- use it!!
						if (iComp < axes[3] - 1) { // Individual components
							if (type == "isotropic")
								array[index] = galaxy.IC_iso_emiss[iComp].d2[iR][iZ].s[iP];
							else
								array[index] = galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP];
						} else { // Total emissivity
							double sum = 0;
							for (unsigned int iC = 0; iC < axes[3]-1; ++iC)
								if (type == "isotropic")
									sum += galaxy.IC_iso_emiss[iC].d2[iR][iZ].s[iP];
								else
									sum += galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP];
							array[index] = sum;
						}
						array[index] *= galaxy.E_gamma[iP]*galaxy.E_gamma[iP];
					} //iR
				} //iZ
			} //iP
		} //iComp
	} //2D

	if(3==galaxy.n_spatial_dimensions) {
		for (unsigned int iComp = 0; iComp < axes[4]; ++iComp) {
			for (unsigned int iP = 0; iP < axes[3]; ++iP) {
				for (unsigned int iZ = 0; iZ < axes[2]; ++iZ) {
					for (unsigned int iY = 0; iY < axes[1]; ++iY) {
						for (unsigned int iX = 0; iX < axes[0]; ++iX) {
							unsigned int index = iComp*axes[3]*axes[2]*axes[1]*axes[0] + iP*axes[2]*axes[1]*axes[0] + iZ*axes[1]*axes[0] + iY*axes[0] + iX; // explicit encoding -- use it!!
							if (iComp < axes[4] - 1) { // Individual components
								if (type == "isotropic")
									array[index] = galaxy.IC_iso_emiss[iComp].d3[iX][iY][iZ].s[iP];
								else
									array[index] = galaxy.IC_aniso_emiss->d3[iX][iY][iZ].s[iP];
							} else { // Total emissivity
								double sum = 0;
								for (unsigned int iC = 0; iC < axes[4]-1; ++iC)
									if (type == "isotropic")
										sum += galaxy.IC_iso_emiss[iC].d3[iX][iY][iZ].s[iP];
									else
										sum += galaxy.IC_aniso_emiss->d3[iX][iY][iZ].s[iP];
								array[index] = sum;
							}
							array[index] *= galaxy.E_gamma[iP]*galaxy.E_gamma[iP];
						} //iX
					} //iY
				} //iZ
			} //iP
		} //iComp
	} //3D


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
	} //2D

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
	} //3D

	fits_close_file(fptr, &status);
	fits_report_error(stderr, status);

	INFO("Exit");

}
