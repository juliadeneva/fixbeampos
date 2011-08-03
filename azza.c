#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

#define CTD (360.0/(M_PI*2.0))
#define CTR (M_PI*2.0/360.0)
#define AOLAT   18.353806         /* degrees north */

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

/* returns the azimuth and za for a source at apparent ra, dec and lst
   ra (hours)
   dec (degrees)
   lst (hours)
*/

int compute_azza( ra, dec, lst, paz, pza  ) 
double ra, dec, lst, *paz, *pza;
{
  double ha, az, el, sindec, cosdec, sinlat, coslat,
                cosha, sinel, cosel, cosazel, cosaz;


/* Calculate hour angle and convert to degrees (use convention
   that hour angle ranges between +/- 12 hours)    */

  ha = lst - ra;
  if( ha > 12.0 )
    ha = ha - 24.0;
  else if( ha < -12.0 )
    ha = ha + 24.0;
  ha = ha * 15.0;

/* Compute trig quantities needed; make angle calculations double
   precision to cut down on round-off-induced overflows near the
   zenith and the meridian.    */

  sindec = sin(dec*CTR);
  cosdec = cos(dec*CTR);
  sinlat = sin(AOLAT*CTR);
  coslat = cos(AOLAT*CTR);
  cosha  = cos(ha*CTR);

/* Compute elevation (trap for blow-ups near the zenith):  */

  sinel = sindec * sinlat + cosdec * cosha * coslat;

  if (sinel >= -1.0 && sinel <= 1.0)
    el = asin(sinel)*CTD;
  else
    printf("elevation too high");

/* Compute azimuth (trap for blow-ups near HA = 0):    */

  cosel = cos(el*CTR);
  cosazel = sindec * coslat - cosdec * cosha * sinlat;
  cosaz = (1.0/cosel) * cosazel;
  if (cosaz >= -1.0 && cosaz <= 1.0)
    az = acos(cosaz)*CTD;
  else {
    if (dec <= AOLAT)
      az = 180.0;
    else
      az = 0.0;
  }

/* Resolve quadrant ambiguity of azimuth:    */

  if (ha > 0.0)
    az = 360.0 - az;

  *pza = 90.0 - el;

  az += 180.0; /* arecibo az is 180 away from azel dish */
  if ( az > 360.0 )
    az -= 360.0;

  *paz = az;
  /* printf("az %f dec %f ha %f\n", az, dec, ha ); */
  return 0;
}

/* 
  thanks to Avinash Deshspande, he computed these mixed up 
  equations
*/

compute_zaha( az, dec, pza, pha )
double az;
double dec;
double *pza;
double *pha;
{
  double sinaz, cosaz, sindec, cosdec, sinlat, coslat;
  double sinelp, sinelm, sinel, cosel, sinha;
  double el, ha, za, sq, zam, zap;

  sinaz  = sin(az*CTR);
  cosaz  = cos(az*CTR);
  sindec = sin(dec*CTR);
  cosdec = cos(dec*CTR);
  sinlat = sin(AOLAT*CTR);
  coslat = cos(AOLAT*CTR);

  sq = cosdec*cosdec - coslat*coslat*sinaz*sinaz;
  if( sq < 0.0 )
    return;
  sinelp = (sinlat*sindec + coslat * cosaz * 
    sqrt( sq ) ) /
    ( 1.0 - coslat*coslat * sinaz*sinaz );

  sinelm = (sinlat*sindec - coslat * cosaz * 
    sqrt( sq ) ) /
    ( 1.0 - coslat*coslat * sinaz*sinaz );

  zam = 90.0 - deg_asin( sinelm );
  zap = 90.0 - deg_asin( sinelp );

/* find a za in the right range */

  if( zam > 0.0 && zam < 20.0 ) {
    sinel = sinelm;
    za = zam;
  } else {
    sinel = sinelp;
    za = zap;
  }
    
  cosel = (sindec - sinel * sinlat ) / coslat / cosaz ;

  sinha =  -1.0 * cosel * sinaz / cosdec;
  ha = deg_asin( sinha );
  ha /= 15;

  *pza = za;
  *pha = ha;
}

/* compute radec offsets from az/za offsets and positions */

int off_hadec( daz, dza, az, za, ha, dec, pdha, pddec )
double daz;
double dza; 
double az; 
double za; 
double ha; 
double dec; 
double *pdha;
double *pddec;
{
  double cosza, coslat, cosaz, cosdec;
  double sinza, sinlat, sinaz;
  double cosha, sinha, ddec, hacd;

   daz = CTR*daz;
   dza = CTR*dza;
   az = CTR*az;
   za = CTR*za;
   ha =  ha*M_PI/12.0;
   dec = CTR*dec;
  
   cosza = cos(za);
   coslat = cos(AOLAT);
   cosaz = cos(az);
   sinza = sin(za);
   sinlat = sin(AOLAT);
   sinaz = sin(az);
   cosdec = cos(dec);
   cosha = cos(ha);
   cosdec = cos(dec);

   ddec = (dza * ( cosza *coslat * cosaz - sinza * sinlat ) - 
          daz * ( coslat * sinaz ) ) / cosdec; 

   *pddec = CTD*ddec;

   hacd = -1.0 * (cosha * cosaz + sinha * sinaz *sinlat ) * daz -
                 (cosha * sinaz * cosza - 
                   sinha * ( sinza * coslat + cosza * sinlat * cosaz )) * dza;

   *pdha = hacd / cosdec;
}


