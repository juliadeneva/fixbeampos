#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "psrfits.h"

//Mock data has correct coordinates starting with v3.43
#define HDRVERGOOD 3.43 //both main header and rows fixed
#define HDRVERFIX 3.429 //only main header fixed, not rows
#define NUMFILES 7

float ddmmss2deg(char* ddmmss, int issign);
char* deg2ddmmss(float deg, int issign);
void alfa_position(double ra, double dec, double lst, double epoch, double angle, double off1, double off2, int beam, double *pra, double *pdec, double *paz, double *pza );
void glgb(double radeg, double decdeg, double* gl, double* gb);

int main(int argc, char** argv)
{
  FILE* outfile;
  int status, ii, beam, beamnum, numrows=0, kk, rowcount;
  struct psrfits pfin;
  char *pc1, *pc2, *ibeam, *beamrastr, *beamdecstr;
  fitsfile *infits, *outfits;
  char tmp[24], hdrver[24], histstr[200];
  double current_epoch;
  double rahh, decdd, hdrverf;
  double beamrahh[NUMFILES], beamdecdd[NUMFILES], beamaz, beamza;
  time_t currtime;
  float* tel_az[NUMFILES], *tel_zen[NUMFILES];
  double *ra_sub[NUMFILES], *dec_sub[NUMFILES], *glon_sub[NUMFILES], *glat_sub[NUMFILES];

  if (argc < 2 || argc > 8) {
    printf("Usage: fixbeampos [beam0 psrfits file] [beam 1 - beam 6 files]\n");
    exit(1);
  }


  // Open the input files
  status = 0;  //fits_close segfaults if this is not initialized
  fprintf(stderr, "Reading input data from:\n");
  for (ii = 0; ii < argc-1; ii++) {
    
    fprintf(stderr,"  '%s'\n", argv[ii+1]);
    //Get the file basename and number from command-line argument
    //(code taken from psrfits2fil)
    pc2 = strrchr(argv[ii+1], '.');      // at .fits
    *pc2 = 0;               // terminate string
    pc1 = pc2 - 1;
    while ((pc1 >= argv[ii+1]) && isdigit(*pc1))
      pc1--;
    if (pc1 <= argv[ii+1]) {     // need at least 1 char before filenum
      puts("Illegal input filename. must have chars before the filenumber");
      exit(1);
    }
    pc1++;                  // we were sitting on "." move to first digit
    pfin.filenum = atoi(pc1);
    pfin.fnamedigits = pc2 - pc1;   // how many digits in filenumbering scheme.
    *pc1 = 0;               // null terminate the basefilename
    strcpy(pfin.basefilename, argv[ii+1]);
    pfin.initialized = 0;   // set to 1 in  psrfits_open()
    pfin.status = 0;
    
    //Get the beam number from the file name
    ibeam = strrchr(argv[ii+1], 'b');
    ibeam = ibeam+1;
    *(ibeam+1) = 0;  //terminate string
    beamnum = atoi(ibeam);
    //fprintf(stderr,"ibeam: %s\n", ibeam);
    
    if (beamnum == 0 && ii != 0) {
      fprintf(stderr, "Stray beam 0 file in list, will skip!\n");
      fprintf(stderr, "First beam 0 file encountered: %s\n", argv[1]);
      fprintf(stderr, "This file: %s\n", argv[ii+1]);
      continue;
    }

    if (ii == 0) {
      //Reading the beam 0 file
      
      if (beamnum != 0) {
	fprintf(stderr, "File is not from ALFA beam 0!\n");
	exit(1);
      }
      
      //Open the existing psrfits file
      if (psrfits_open(&pfin, READONLY) != 0) {
	    fprintf(stderr, "error opening file\n");
	    fits_report_error(stderr, pfin.status);
	    exit(1);
      }
      infits = pfin.fptr;

      //Allocate various arrays needed by psrfits_read_subint
      pfin.sub.dat_freqs = (float *) malloc(sizeof(float) * pfin.hdr.nchan);
      pfin.sub.dat_weights = (float *) malloc(sizeof(float) * pfin.hdr.nchan);
      pfin.sub.dat_offsets =
        (float *) malloc(sizeof(float) * pfin.hdr.nchan * pfin.hdr.npol);
      pfin.sub.dat_scales =
        (float *) malloc(sizeof(float) * pfin.hdr.nchan * pfin.hdr.npol);

      //Move to main HDU
      fits_movabs_hdu(infits, 1, NULL, &status);

      //Check that receiver is ALFA and backend is Mock
      if(strcmp(pfin.hdr.frontend, "alfa") != 0) {
	    fprintf(stderr,"File not from frontend ALFA!\n");
	    exit(1);
      }
      if(strcmp(pfin.hdr.backend, "pdev") != 0) {
	    fprintf(stderr,"File not from backend Mock!\n");
	    exit(1);
      }

      // Find the epoch in years after epoch 2000.0
      current_epoch = 2000.0 + (pfin.hdr.MJD_epoch - 51544.0) / 365.25;
      
      //Find RA & DEC as decimal hours and degrees, respectively
      //RA is a string of the form HH:MM:SS.SSSS
      //DEC is a string of the form +DD:MM:SS.SSSS
      //rahh = ddmmss2deg(rastr,0);
      //printf("rahh: %f\n", rahh);
      //decdd = ddmmss2deg(decstr,1);
      //printf("decdd: %f\n", decdd);

      //Move to SUBINT table
      fits_movnam_hdu(infits, BINARY_TBL, "SUBINT", 0, &pfin.status);
      //Get the number of rows
      numrows = pfin.rows_per_file;
      fprintf(stderr,"numrows: %d\n",numrows);
      //Allocate coordinate arrays for all beams
      for(kk=0; kk<NUMFILES; kk++) {
        ra_sub[kk] = (double*) malloc(sizeof(double) * numrows);
        dec_sub[kk] = (double*) malloc(sizeof(double) * numrows);
        glon_sub[kk] = (double*) malloc(sizeof(double) * numrows);
        glat_sub[kk] = (double*) malloc(sizeof(double) * numrows);
        tel_az[kk] = (float*) malloc(sizeof(float) * numrows);
        tel_zen[kk] = (float*) malloc(sizeof(float) * numrows);
      }
	
      
      outfile = fopen("beampos.out","w");
      rowcount = 0;
      while (psrfits_read_subint(&pfin) == 0) {
        fprintf(stderr, "Working on row %d\r", rowcount+1);

      if(rowcount == 0)
        fprintf(outfile,"# b0az b0za b0ra b0dec b1az b1za b1ra b1dec, etc\n");

	//printf("row: %d  feed_ang: %f  tel_az: %f tel_zen: %f ra_sub: %f dec_sub: %f\n", rowcount, pfin.sub.feed_ang, pfin.sub.tel_az, pfin.sub.tel_zen, pfin.sub.ra, pfin.sub.dec);
	//NOTE: ra_sub is in degrees, while alfa_position wants the ra in hours
	
	// Save beam 0 position for debugging
	//Save this row's RA, DEC, AZ & ZA
	tel_az[0][rowcount] = pfin.sub.tel_az;
	tel_zen[0][rowcount] = pfin.sub.tel_zen;
	ra_sub[0][rowcount] = pfin.sub.ra;
	dec_sub[0][rowcount] = pfin.sub.dec;
	glon_sub[0][rowcount] = pfin.sub.glon;
	glat_sub[0][rowcount] = pfin.sub.glat;

	decdd = pfin.sub.dec;
	rahh = pfin.sub.ra/360.0 * 24.0;

	for(beam=1; beam<NUMFILES; beam++) {
        alfa_position(rahh, decdd, pfin.sub.lst/3600.0, 
                      current_epoch, pfin.sub.feed_ang, 
                      0.0, 0.0, beam, &beamrahh[beam], &beamdecdd[beam], 
                      &beamaz, &beamza);
	  
	  //Save this row's RA, DEC, AZ & ZA, GL & GB
	  tel_az[beam][rowcount] = (float)beamaz;
	  tel_zen[beam][rowcount] = (float)beamza;
	  ra_sub[beam][rowcount] = beamrahh[beam]/24.0 * 360.0;
	  dec_sub[beam][rowcount] = beamdecdd[beam];

	  glgb(ra_sub[beam][rowcount], dec_sub[beam][rowcount], 
           &glon_sub[beam][rowcount],&glat_sub[beam][rowcount]);

	  //Convert beam RA & DEC from 1st row to strings to put in side 
	  //beam main HDU
	  if(rowcount == 0) {
	    beamrastr = deg2ddmmss(beamrahh[beam],0);
	    beamdecstr = deg2ddmmss(beamdecdd[beam],1);
	  }
	}
	      
	//Print positions for all beams to debugging file
	for(beam=0; beam<NUMFILES; beam++) {
        fprintf(outfile,"%f %f %f %f %f %f ",tel_az[beam][rowcount],tel_zen[beam][rowcount],ra_sub[beam][rowcount],dec_sub[beam][rowcount],glon_sub[beam][rowcount],glat_sub[beam][rowcount]);
	}
	fprintf(outfile,"\n");
    
	rowcount++;
      }

      //Print fits file name & position for Patrick's script
      fprintf(stdout,"%s  %s  Old  %s  %s  %f  %f\n",pfin.filename,ibeam,pfin.hdr.ra_str, pfin.hdr.dec_str, ra_sub[beamnum][0],dec_sub[beamnum][0]);
      fprintf(stdout,"%s  %s  New  %s  %s  %f  %f\n",pfin.filename,ibeam,pfin.hdr.ra_str, pfin.hdr.dec_str, ra_sub[beamnum][0],dec_sub[beamnum][0]);

      fclose(outfile);
      fits_close_file(infits, &pfin.status);

    } else {

      //Open side beam file for writing
      if (psrfits_open(&pfin, READWRITE) != 0) {
	    fprintf(stderr, "Error opening file\n");
        //fits_report_error(stderr, pfin.status); //this exits
        continue;
      }
      outfits = pfin.fptr;
      
      //Go to main HDU
      fits_movabs_hdu(outfits, 1, NULL, &status);
      
      //Check version number
      fits_read_key(outfits, TSTRING, "HDRVER", hdrver, NULL, &pfin.status);
      hdrverf = atof(hdrver);
      //printf("hdrver: %f\n",hdrverf);
      
      //Convert this beam's RA & DEC to strings
      beamrastr = deg2ddmmss(beamrahh[beamnum],0);
      beamdecstr = deg2ddmmss(beamdecdd[beamnum],1);
      //printf("beam: %d  ra: %s  dec: %s\n", ii,beamrastr,beamdecstr);
      

      //UNCOMMENT FOR RELEASE!!!
      if (fabs(hdrverf-HDRVERGOOD) < 0.001) {
        printf("HDRVER = %5.3f, file does not need fixing.\n",hdrverf);
        continue;
      }

      //Replace RA & DEC with correct values
      fits_update_key(outfits, TSTRING, "RA", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "DEC", beamdecstr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STT_CRD1", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STT_CRD2", beamdecstr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STP_CRD1", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STP_CRD2", beamdecstr, NULL, &pfin.status);
      
      //Update version number
      sprintf(hdrver,"%5.3f",HDRVERGOOD);
      fits_update_key(outfits, TSTRING, "HDRVER", hdrver, NULL, &pfin.status);
      //Add field for fix date
      time(&currtime);                                                     
      strftime(tmp, sizeof(tmp)-1, "%Y-%m-%d",localtime(&currtime));
      //printf("FIXDATE: %s \n",tmp);
      // Add a HISTORY line with the date fixed and the old coords
      sprintf(histstr, "Positional information changed by fixbeampos on %s.", tmp);
      fits_write_history(outfits, histstr, &pfin.status);
      sprintf(histstr, "Original RA and DEC were '%s', '%s'",
              pfin.hdr.ra_str, pfin.hdr.dec_str);
      fits_write_history(outfits, histstr, &pfin.status);
      sprintf(histstr, "Git hash of fixbeampos: %s", GITHASH);
      fits_write_history(outfits, histstr, &pfin.status);

      //Move to start of SUBINT table
      fits_movnam_hdu(outfits, BINARY_TBL, "SUBINT", 0, &pfin.status);

      //Write corrected positions into rows
      for(rowcount=1; rowcount<=numrows; rowcount++) {
        fprintf(stderr, "Correcting row %d\r", rowcount);
        fits_write_col(outfits, TDOUBLE, pfin.subcols.ra_sub, rowcount, 1, 1, &ra_sub[beamnum][rowcount-1], &pfin.status);
        fits_write_col(outfits, TDOUBLE, pfin.subcols.dec_sub, rowcount, 1, 1, &dec_sub[beamnum][rowcount-1], &pfin.status);
        fits_write_col(outfits, TDOUBLE, pfin.subcols.glon_sub, rowcount, 1, 1, &glon_sub[beamnum][rowcount-1], &pfin.status);
        fits_write_col(outfits, TDOUBLE, pfin.subcols.glat_sub, rowcount, 1, 1, &glat_sub[beamnum][rowcount-1], &pfin.status);
        fits_write_col(outfits, TFLOAT, pfin.subcols.tel_az, rowcount, 1, 1, &tel_az[beamnum][rowcount-1], &pfin.status);
        fits_write_col(outfits, TFLOAT, pfin.subcols.tel_zen, rowcount, 1, 1, &tel_zen[beamnum][rowcount-1], &pfin.status);
      }

      /*
      //Open debugging file
      sprintf(outfilename,"%s%d%s","beampos",ii,".out.after");
      outfile = fopen(outfilename,"w");

      //Move to start of SUBINT table
      fits_movnam_hdu(outfits, BINARY_TBL, "SUBINT", 0, &pfin.status);

      //Print out positions in all rows for debugging
      rowcount = 0;
      while (psrfits_read_subint(&pfin) == 0) {
	fprintf(stderr, "Reading row %d\r", rowcount+1);

	fprintf(outfile,"%f %f %f %f %f %f\n",pfin.sub.tel_az,pfin.sub.tel_zen,pfin.sub.ra,pfin.sub.dec,pfin.sub.glon,pfin.sub.glat);

	rowcount++;
      }
      fclose(outfile);		     
      */

      //Print fits file name & position for Patrick's script
      fprintf(stdout,"%s  %s  Old  %s  %s  %f  %f\n",pfin.filename, ibeam, pfin.hdr.ra_str, pfin.hdr.dec_str, 15.0*ddmmss2deg(pfin.hdr.ra_str,0), ddmmss2deg(pfin.hdr.dec_str,1));
      fprintf(stdout,"%s  %s  New  %s  %s  %f  %f\n",pfin.filename, ibeam, beamrastr, beamdecstr, 15.0*beamrahh[beamnum], beamdecdd[beamnum]);

      
      fits_close_file(outfits, &pfin.status);
    }//end if(ii == 0) else

  }//end for (ii = 0; ii < NUMFILES; ii++)

  free(pfin.sub.dat_freqs);
  free(pfin.sub.dat_weights);
  free(pfin.sub.dat_offsets);
  free(pfin.sub.dat_scales);

  for(kk=0; kk<NUMFILES; kk++) {
    free(ra_sub[kk]);
    free(dec_sub[kk]);
    free(glon_sub[kk]);
    free(glat_sub[kk]);
    free(tel_az[kk]);
    free(tel_zen[kk]);
  }
  
  return 0;

}


//Convert from sexagesimal to decimal degrees (or decimal hours for RA).
//The 2nd argument is nonzero if the first character in the string contains
//a sign (e.g. for DEC)
float ddmmss2deg(char* ddmmss, int issign)
{
  int dd, mm;
  float ss, deg;
  char ddmmsstemp[16];
  
  strcpy(ddmmsstemp, ddmmss);

  ddmmsstemp[2+issign] = '\0';
  dd = atoi(ddmmsstemp);
  ddmmsstemp[5+issign] = '\0';
  mm = atoi(&ddmmsstemp[3+issign]);
  ss = atof(&ddmmsstemp[6+issign]);
  deg = (float)dd + (float)mm/60.0 + ss/3600.0;

  if(issign && ddmmsstemp[0] == '-')
    deg = -deg;
  
  return deg;

}

//Convert from decimal to sexagesimal. The 2nd argument determines whether
//the sign will be included in the output string.
char* deg2ddmmss(float deg, int issign)
{
  int dd, mm, sign;
  float ss, tmp;

  char* ddmmss = (char *) malloc(sizeof(char)*24);

  if(deg > 0.0)
    sign = 1;
  else
    sign = -1;

  deg = fabs(deg);

  dd = floor(deg);
  tmp = (deg - (float)dd) * 60.0;
  mm = floor(tmp);
  ss = (deg - (float)dd - (float)mm/60.0)*3600.0;

  if (issign) {
    if(sign == 1)
      sprintf(ddmmss, "+%.2d:%.2d:%07.4f",dd,mm,ss);
    else
      sprintf(ddmmss, "-%.2d:%.2d:%07.4f",dd,mm,ss);
  } 
  else
    sprintf(ddmmss, "%.2d:%.2d:%07.4f",dd,mm,ss);

  return ddmmss;
  
}
