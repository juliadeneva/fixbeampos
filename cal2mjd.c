#include <stdio.h> 
#include <math.h> 

/*---------------------------------------------------------------------*
 * return modified julian date given a gregorian calendar - based on slalib
 * routine sla_cldj which itself is based on Hatcher (1984) QJRAS 25, 53-55 
 * creation date 1999/07/10 [mjd=51639] (dunc@naic.edu) 
 *---------------------------------------------------------------------*/
int cal2mjd(int iy, int im, int id) 
{
  int leap;
  double tmp;

  //fprintf(stderr, "cal2mjd: year-month-day: %d-%d-%d\n",iy,im,id);

  /* month lengths in days for normal and leap years */
  static int mtab[2][13] = {
    {0,31,28,31,30,31,30,31,31,30,31,30,31},
    {0,31,29,31,30,31,30,31,31,30,31,30,31}
  };

  /*validate year*/
  if (iy<-4699) {
    fprintf(stderr,"Invalid year passed to cal2mjd\n");
    return(0.0);
  } else {
    /* validate month */
    if (im<1 || im>12) {
      fprintf(stderr,"Invalid month passed to cal2mjd\n");
      return(0.0);
    } else {
      /* allow for leap year */
      leap = iy%4 == 0 && iy%100 != 00 || iy%400 == 0;
      /* validate day */
      if (id<1 || id>mtab[leap][im]) {
	fprintf(stderr,"Invalid day passed to cal2mjd\n");
        return(0.0);
      }
    }
  }
  
  return (1461*(iy-(12-im)/10+4712))/4 + (5+306*((im+9)%12))/10 - (3*((iy-(12-im)/10+4900)/100))/4 + id - 2399904;

}

