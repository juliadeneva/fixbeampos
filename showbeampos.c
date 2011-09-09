#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
//#include <time.h>

#define NUMFILES 7

int main(int argc, char** argv)
{
  FILE* infile;
  char fitsname[128];
  float radeg, decdd, epoch, lsthh, feed_ang, rahh;
  int beam;
  double beamaz, beamza, beamrahh, beamdecdd;
  char* beamnum, beamnumstr[2];

  if(argc != 2) {
    fprintf(stderr, "Usage: showbeampos [beam0 position list file]\n");
    exit(1);
  }

  infile = fopen(argv[1],"r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s.\n",argv[1]);
    exit(1);
  }

  while(!feof(infile)) {
    fscanf(infile,"%s %f %f %f %f %f",fitsname, &radeg, &decdd, &epoch, &lsthh, &feed_ang);
    printf("%-50s %12.4f %8.4f %12.4f %12.6f %8.4f\n",fitsname, radeg, decdd, epoch, lsthh, feed_ang);

    //Find beam number in file name
    beamnum = strstr(fitsname, ".b");
    
    for(beam=1; beam<NUMFILES; beam++){

      //Construct side beam file name
      sprintf(beamnumstr, "%d",beam);
      *(beamnum+2) = beamnumstr[0];
      //printf("%s \n",fitsname);

      //Calculate side beam position
      rahh = radeg/360.0 * 24.0;      

      alfa_position((double)rahh, (double)decdd, (double)lsthh, (double)epoch, (double)feed_ang, 0.0, 0.0, beam, &beamrahh, &beamdecdd, &beamaz, &beamza);

      printf("%-50s %12.4f %8.4f %12.4f %12.6f %8.4f\n",fitsname, beamrahh/24.0*360.0, beamdecdd, epoch, lsthh, feed_ang);
    }
    
    
  }

  return 0;
}
