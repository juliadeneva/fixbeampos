/*
Convert from RA & DEC to Gl & Gb

*/
#include <stdio.h>
#include <math.h>
#define PI		3.14159265358979323846
#define NGPDEC 27.118 //DEC of North Galactic pole
//#define NGPDEC 27.128 //DEC of North Galactic pole J2000 (K.Lang, 'Astrophysical Formulae, vol. 2')
//#define NGPRA 192.859481
#define A0DEG 282.86
#define L0DEG 32.93

//Arguments are in deg
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
//double  sin(), cos(), tan();

void glgb(double radeg, double decdeg, double* gl, double* gb)
{

  double sinb, cosl_l0, sinl_l0;

  sinb = deg_sin(decdeg)*deg_sin(NGPDEC) - deg_cos(decdeg)*deg_cos(NGPDEC)*deg_sin(radeg-A0DEG);

  *gb = deg_asin(sinb);

  cosl_l0 = deg_cos(radeg-A0DEG)*deg_cos(decdeg)/deg_cos(*gb);
  sinl_l0 = (deg_sin(decdeg)*deg_cos(NGPDEC)+deg_cos(decdeg)*deg_sin(NGPDEC)*deg_sin(radeg-A0DEG))/deg_cos(*gb);


  if(cosl_l0 > 0.0)
      {
	//fprintf(stderr, "Cos l-l0 > 0\n");
	if(sinl_l0 > 0.0)
	  {
	    //fprintf(stderr, "Sin l-l0 > 0\n");
	    *gl = deg_asin(sinl_l0) + L0DEG;
	  }
	else
	  {
	    //fprintf(stderr, "Sin l-l0 < 0 \n");
	    *gl = 360.0+deg_asin(sinl_l0) + L0DEG;
	  }
      }
    else
      {
	//fprintf(stderr, "Cos l-l0 < 0 \n");

	if(sinl_l0 > 0.0)
	  {
	    //fprintf(stderr, "Sin l-l0 > 0\n");
	    *gl = 180.0-deg_asin(sinl_l0) + L0DEG;
	  }
	else
	  {
	    //fprintf(stderr, "Sin l-l0 < 0\n");
	    *gl = 360.0+deg_asin(sinl_l0) + L0DEG;
	  }
      }

  
}
