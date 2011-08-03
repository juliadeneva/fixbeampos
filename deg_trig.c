/*
 * Cactus File @(#)deg_trig.c	1.1
 *         SID 1.1
 *        Date 02/04/98
 */

static char SccsId[] = "@(#)deg_trig.c	1.1\t02/04/98";

#include <math.h>

#define PI		3.14159265358979323846
#define DEG_TO_RAD	PI/180.0
#define RAD_TO_DEG	180.0/PI

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

/* Calculates trigonometric fn's with 'angle' in degrees, not radians.  */

double deg_sin(angle)
double angle;
{
  double sin();

  return( sin(DEG_TO_RAD * angle) );
}

double deg_cos(angle)
double angle;
{
  double cos();

  return( cos(DEG_TO_RAD * angle) );
}

double deg_tan(angle)
double  angle;
{
  double tan();

  return( tan(DEG_TO_RAD * angle) );
}

double deg_asin(value)
double  value;
{
  double asin();

  return( RAD_TO_DEG*asin(value) );
}

double deg_acos(value)
double  value;
{
  double acos();

  return( RAD_TO_DEG*acos(value) );
}

double deg_atan(value)
double  value;
{
  double atan();

  return( RAD_TO_DEG*atan(value) );
}


