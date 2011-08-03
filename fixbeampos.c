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

int obs2mjd(char* dt);
float ddmmss2deg(char* ddmmss, int issign);
char* deg2ddmmss(float deg, int issign);
void alfa_position(double ra, double dec, double lst, double epoch, double angle, double off1, double off2, int beam, double *pra, double *pdec );

int main(int argc, char** argv)
{
  FILE* outfile;
  int status, ii, beam, beamnum;
  struct psrfits pfin;
  char *pc1, *pc2, *ibeam, *beamrastr, *beamdecstr;
  fitsfile *infits, *outfits;
  char rastr[24], decstr[24], obs_date[24], datecat[24], tmp[24], hdrver[24];
  char frontend[24], backend[24];
  int stt_imjd, stt_smjd;
  float stt_offs, stt_lst, start_mjd, start_lst, current_time, mjd, secs;
  float rahh, decdd, hdrverf;
  double beamrahh[7], beamdecdd[7];
  time_t currtime;

  if (argc != 2) {
    printf("Usage: fixbeampos [beam0 psrfits file]\n");
    printf("(All psrfits files should be in the current directory.)\n");
    exit(1);
  }

  // Open the input files
  status = 0;  //fits_close segfaults if this is not initialized
  printf("Reading input data from:\n");
  for (ii = 0; ii < NUMFILES; ii++) {
    if (ii == 0) {
      printf("  '%s'\n", argv[ii+1]);
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
      
      //Reading the beam 0 file
      
      if (beamnum != 0) {
	printf("File is not from ALFA beam 0!\n");
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
      fits_read_key(infits, TSTRING, "FRONTEND", frontend, NULL, &pfin.status);
      if(strcmp(frontend, "alfa") != 0) {
	fprintf(stderr,"File not from frontend ALFA!\n");
	exit(1);
      }
      fits_read_key(infits, TSTRING, "BACKEND", backend, NULL, &pfin.status);
      if(strcmp(backend, "pdev") != 0) {
	fprintf(stderr,"File not from backend Mock!\n");
	exit(1);
      }

      //Read various fields needed for calculation of the side beam positions
      fits_read_key(infits, TSTRING, "DATE-OBS", obs_date, NULL, &pfin.status);
      //printf("obs_date: %s\n", obs_date);
      fits_read_key(infits, TSTRING, "RA", rastr, NULL, &pfin.status);
      //printf("rastr: %s\n",rastr);
      fits_read_key(infits, TSTRING, "DEC", decstr, NULL, &pfin.status);
      //printf("decstr: %s\n",decstr);
      fits_read_key(infits, TINT, "STT_IMJD", &stt_imjd, NULL, &pfin.status);
      fits_read_key(infits, TINT, "STT_SMJD", &stt_smjd, NULL, &pfin.status);
      fits_read_key(infits, TFLOAT, "STT_OFFS", &stt_offs, NULL, &pfin.status);
      //printf("status: %d\n", pfin.status);
      fits_read_key(infits, TFLOAT, "STT_LST", &stt_lst, NULL, &pfin.status);
      //printf("stt_imjd: %d  stt_smjd: %d  stt_offs: %f  stt_lst: %f\n",stt_imjd, stt_smjd, stt_offs, stt_lst);
      
      //Find the MJD, LST, and UT days past 2000-01-01 (epoch 2000.0) at start
      secs = (float)stt_smjd + stt_offs; 
      start_mjd = (float)stt_imjd + secs/(24.0*3600.0);
      start_lst = stt_lst/3600.0; 
      
      //Start copy from alfasplit.c
      strcpy(tmp, obs_date);
      tmp[4] = '\0';
      tmp[7] = '\0';
      tmp[10] = '\0';
      strcat(datecat,tmp);
      strcat(datecat,&tmp[5]);
      strcat(datecat,&tmp[8]);
      //printf("datecat: %s\n", datecat);
      mjd = obs2mjd(datecat);
      mjd = mjd - obs2mjd("20000101");
      current_time = 2000.0 + (mjd + secs/(3600.0*24.0))/365.25; //ut
      //End copy from alfasplit.c
      
      //Find RA & DEC as decimal hours and degrees, respectively
      //RA is a string of the form HH:MM:SS.SSSS
      //DEC is a string of the form +DD:MM:SS.SSSS
      rahh = ddmmss2deg(rastr,0);
      //printf("rahh: %f\n", rahh);
      decdd = ddmmss2deg(decstr,1);
      //printf("decdd: %f\n", decdd);

      //Move to SUBINT table
      fits_movnam_hdu(infits, BINARY_TBL, "SUBINT", 0, &pfin.status);
      //printf("status: %d\n", pfin.status);
      //Read ALFA rotation angle from first row
      psrfits_read_subint(&pfin, 1);
      //printf("feed_ang: %f pos_ang: %f  par_ang: %f\n", pfin.sub.feed_ang, pfin.sub.pos_ang, pfin.sub.par_ang);
      
      outfile = fopen("beampos.out","w");
      //Calculate the side beam positions
      for(beam=1; beam<NUMFILES; beam++) {
	alfa_position((double)rahh, (double)decdd, (double)start_lst, (double)current_time, (double)pfin.sub.feed_ang, 0.0, 0.0, beam, &beamrahh[beam], &beamdecdd[beam]);
	
	beamrastr = deg2ddmmss(beamrahh[beam],0);
	beamdecstr = deg2ddmmss(beamdecdd[beam],1);
	
	fprintf(outfile,"beam: %d  rahh: %f  %s  decdd: %f  %s\n", beam,beamrahh[beam], beamrastr, beamdecdd[beam], beamdecstr);
	fprintf(stderr,"beam: %d  rahh: %f  %s  decdd: %f  %s\n", beam,beamrahh[beam], beamrastr, beamdecdd[beam], beamdecstr);
      }
      fclose(outfile);

      fits_close_file(infits, &pfin.status);

    } else {

      //Construct file name for side beam
      ibeam = strrchr(pfin.basefilename, 'b');
      pc2 = ibeam+2; //portion of basefilename after beam number
      ibeam = ibeam+1;
      *(ibeam) = 0;  //terminate string
      beamnum = ii;
      sprintf(pfin.basefilename,"%s%d%s",pfin.basefilename,ii,pc2);
      //printf("beam: %d pfin.basefilename: %s\n",ii,pfin.basefilename);

      pfin.initialized = 0;   // set to 1 in  psrfits_open()
      pfin.status = 0;

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
      
      //UNCOMMENT FOR RELEASE!!!
      if(fabs(hdrverf-HDRVERFIX) < 0.001) {
	printf("HDRVER = %5.3f, file does not need fixing.\n",hdrverf);
	continue;
      }
      
      //Convert this beam's RA & DEC to strings
      beamrastr = deg2ddmmss(beamrahh[ii],0);
      beamdecstr = deg2ddmmss(beamdecdd[ii],1);
      //printf("beam: %d  ra: %s  dec: %s\n", ii,beamrastr,beamdecstr);
      
      //Replace RA & DEC with correct values
      fits_update_key(outfits, TSTRING, "RA", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "DEC", beamdecstr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STT_CRD1", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STT_CRD2", beamdecstr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STP_CRD1", beamrastr, NULL, &pfin.status);
      fits_update_key(outfits, TSTRING, "STP_CRD2", beamdecstr, NULL, &pfin.status);
      
      //Update version number
      sprintf(hdrver,"%5.3f",HDRVERFIX);
      fits_update_key(outfits, TSTRING, "HDRVER", hdrver, NULL, &pfin.status);
      //Add field for fix date
      time(&currtime);                                                     
      strftime(tmp,sizeof(tmp)-1,"%Y-%m-%d",localtime(&currtime)); 
      //printf("FIXDATE: %s \n",tmp);
      fits_update_key(outfits, TSTRING, "FIXDATE",tmp,"Side ALFA beam position fix date (YYYY-MM-DD)",&pfin.status);
      
      fits_close_file(outfits, &pfin.status);
    }//end if(ii == 0) else

  }//end for (ii = 0; ii < NUMFILES; ii++)

  free(pfin.sub.dat_freqs);
  free(pfin.sub.dat_weights);
  free(pfin.sub.dat_offsets);
  free(pfin.sub.dat_scales);

  return 0;

}


//Convert from sexagesimal to decimal degrees (or decimal hours for RA).
//The 2nd argument is nonzero if the first character in the string contains
//a sign (e.g. for DEC)
float ddmmss2deg(char* ddmmss, int issign)
{
  int dd, mm;
  float ss, deg;

  ddmmss[2+issign] = '\0';
  dd = atoi(ddmmss);
  ddmmss[5+issign] = '\0';
  mm = atoi(&ddmmss[3+issign]);
  ss = atof(&ddmmss[6+issign]);
  deg = (float)dd + (float)mm/60.0 + ss/3600.0;

  if(issign && ddmmss[0] == '-')
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
