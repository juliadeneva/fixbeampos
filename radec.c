#include <stdlib.h>
#include <math.h>
/*
   compute apparent ra, dec from az, el, lst
   Zombeck 1982 p. 71

   az, el, ra, dec in degrees
   lst in hours

*/

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

#define AOLAT   18.353806         /* degrees north */

/* az za in degrees lst in hours */

compute_radec( az, za, lst, pra, pdec )
double az, za, lst;
double *pra, *pdec;
{
  double sindec, sinaz, sinel, cosaz, cosel;
  double cosdec, dec, ha, cosha, sinlat, coslat;
  double el;

  az += 180.0;
  if( az > 360.0 )
   az -= 360.0;
  el = 90.0 - za;
  sinel = deg_sin(el);
  cosel = deg_cos(el);
  sinaz = deg_sin(az);
  cosaz = deg_cos(az);
  sinlat = deg_sin(AOLAT);
  coslat = deg_cos(AOLAT);

  sindec = sinel * sinlat + cosel * coslat * cosaz;

  if( sindec >= -1.0 && sindec <= 1.0 )
    dec = deg_asin( sindec );
  else
    dec = 90.0;

  cosdec = deg_cos( dec );
  cosha = ( sinel * coslat - cosel * cosaz * sinlat)/cosdec ;

  if( cosha >= -1.0 && cosha <= 1.0 )
    ha = deg_acos( cosha )/15.0;
  else
    ha = 0.0;

/* Resolve quadrant ambiguity */
  
  if( az < 180.0 )
    ha = -ha;

  *pra =  fmod(lst - ha, 24.0);
  if( *pra < 0.0 )
    *pra += 24.0;

  *pdec = dec;
}

