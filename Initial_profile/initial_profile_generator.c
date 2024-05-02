#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define	NDIM 2
#define Sqr(x) ((x)*(x))
///#include <complex.h>
//#include <fftw3.h>
#include"global.h"


void vortex(double v[], double u[], double dy, double dx){

double Uxp, Uxn, Uyp, Uyn;
 int nx, ny;
FILE *outfp4;
char outfile4[25];


FILE *fpvort;
fpvort = fopen("QUIVER_time_0","w");

    for(ny = 0 ; ny < ngrids_Y ; ny ++){

    //sprintf(outfile4,"Vort_profiles_%d", ny+1);
    //outfp4 = fopen(outfile4, "w");

      for(nx = 0 ; nx < ngrids_X ; nx ++){
	if(nx == 0){
	  Uyn = u[(ny)*ngrids_X + 1];
	  Uyp = u[(ny)*ngrids_X + ngrids_X-1];
	}else if(nx == ngrids_X-1){
	  Uyn = u[(ny)*ngrids_X];
	  Uyp = u[(ny)*ngrids_X + ngrids_X-2];
	}else{
	  Uyn = u[(ny)*ngrids_X + nx+1];
	  Uyp = u[(ny)*ngrids_X + nx-1];
	}
	if(ny == 0){
	  Uxn = v[ngrids_X + nx];
	  Uxp = v[(ngrids_Y-1)*ngrids_X + nx];
	}else if(ny == ngrids_Y-1){
	  Uxn = v[nx];
	  Uxp = v[(ngrids_Y-2)*ngrids_X + nx];
	}else{
	  Uxn = v[(ny+1)*ngrids_X + nx];
	  Uxp = v[(ny-1)*ngrids_X + nx];
	}
	n = (ny)*ngrids_X + nx;
	x = (nx-0.5)*dx - halfx;
	y = (ny-0.5)*dy - halfy;
	vorticity[n] = (Uyn - Uyp)/(2.0*dx) - (Uxn - Uxp)/(2.0*dy);
	//if(vorticity[n] <= 0.0)
	//vorticity[n] = 0.0;
	//fprintf(outfp4,"%lf\t %lf\t %lf\t %lf\t %lf %d %d\n",x ,y, v[n], u[n], vorticity[n],ny,nx);
        //fprintf(fpvort,"%lf\t %lf\t %lf\t %lf\t %lf\t %d\n",x ,y, v[n], u[n], omega[n],n);
        fprintf(fpvort,"%lf\t %lf\t %lf\t %lf\t %lf\t %d\n",x ,y, v[n], u[n], vorticity[n],n);
      }
//fclose(outfp4);
    }

fclose(fpvort);
}



int main(){

///fftw_complex  ii;
//ii = 0.0+1.0*I;

twopi = 4.0*asin(1);
int m, N = 100; 
ngrids_Y = 512;
ngrids_X = 512;
L = 10; 

char dummy[25];
double dx, dy; 
halfx = 0.5*L; 
halfy = 0.5*L;
dx = (double)L/(double)ngrids_X;
dy = (double)L/(double)ngrids_Y;

/*
FILE *fp3;
fp3 = fopen("region","r");
fscanf(fp3, "%s %lf %s %lf", dummy, &region[1], dummy, &region[2]);
fclose(fp3);

regionH[1] = 0.5*region[1];
regionH[2] = 0.5*region[2];
*/

double  *psi, *psi1;

//Arrays in real space
omega = (double *)malloc((ngrids_Y*ngrids_X)*sizeof(double));
psi = (double *)malloc((ngrids_Y*ngrids_X)*sizeof(double));
psi1 = (double *)malloc((ngrids_Y*ngrids_X)*sizeof(double));
kx = (double*) malloc(sizeof(double) * ngrids_X);
ky = (double*) malloc(sizeof(double) * ngrids_Y);
v = (double *)malloc((ngrids_Y*ngrids_X)*sizeof(double));
u = (double *)malloc((ngrids_Y*ngrids_X)*sizeof(double));
u1 = (double*)malloc(ngrids_Y * sizeof(double));
in = (double*)malloc(ngrids_X*ngrids_Y * sizeof(double));
vorticity = (double*)malloc(ngrids_X*ngrids_Y * sizeof(double));

/*
//Arrays in fourier space
fomega = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
fpsi = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
ffirst = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
fv = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
fu = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
v1 = (fftw_complex*) fftw_malloc(ngrids_Y * sizeof(fftw_complex));
out = (fftw_complex*) fftw_malloc(ngrids_Y*ngrids_X * sizeof(fftw_complex));
*/

for (ny=0;ny<ngrids_Y;ny++){
  for (nx=0;nx<ngrids_X;nx++){
  n = (ny)*ngrids_X + nx;
 
  psi[n] = 0.0;
  psi1[n] = 0.0;
  }
}




////////////////////////////////
//Generating the values for omega

FILE *fpinitial;
fpinitial = fopen("Initial_0","w");

x= 0.0;
for (i=0;i<ngrids_X;i++){
 x += dx;
 y = 0.0;
 for (j=0;j<ngrids_Y;j++){
 y += dy;
 m = (j)*ngrids_X + i; 
double x1 = (i-0.5)*dx - halfx;
double y1 = (j-0.5)*dy - halfy;
 //if((Sqr(x-halfx) + Sqr(y-halfy))<=1.0){
 if(sqrt(Sqr(x1) + Sqr(y1 + 1.2))<=0.5){
 printf("%lf %lf\n", x1, y1);
 omega[m] = 0.5;
 fprintf(fpinitial,"%lf\t %lf\t %lf\t %lf\t %lf\n",x1,y1,0.0,0.0,omega[m]);
 }else if(sqrt(Sqr(x1) + Sqr(y1 - 1.2))<=0.5){
 omega[m] = 0.5;
 }else{
 omega[m] = 0.0;
 fprintf(fpinitial,"%lf\t %lf\t %lf\t %lf\t %lf\n",x1,y1,0.0,0.0,omega[m]);
 }
 }
 //fprintf(fpinitial,"\n");
}

fclose(fpinitial);

/*

//////////////////////////////////
//////////////////////////////////


getarray(omega);
for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){
 r = j*ngrids_X+i;
 fomega[r] = out[r];

 }
}
/////////////////////////////////
/////////////////////////////////

////////////////////////////////
//Initiasation of psi

for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){
 r = j*ngrids_X+i;
 psi[r] = 0.0;
 }
}


//////////////////////////////
//////////////////////////////
//fft of psi


getarray(psi);
for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){
 r = j*ngrids_X+i;
 fpsi[r] = out[r];

 }
}

*/
///////////////////////////
///////////////////////////



//wavevector defination

for(i=0;i<ngrids_Y;i++){
 ky[i] = (twopi/Sqr(L))*(i);
  for(j=0;j<0.5*(ngrids_X);j++){
  kx[j] = (twopi/(Sqr(L)))*(j);
  kx[ngrids_X-j-1] = -(twopi/(Sqr(L)))*(j+1);
  //printf("%lf\t %lf\t %d %d\n",kx[j],kx[ngrids_X-j-1],j,ngrids_X-j-1);
 }
}


//////////////////////////////
//////////////////////////////
int tend = 1, t_end1 = 1000, t;
double first, second;
double sorcoeff = 1.0;


for(t=1;t<=tend;t++){


int iter = 0;
double error = 1.0;
double tol = 1.0e-6;
//while(error>tol && iter < t_end1)
while(error>tol)
{
double sumerr = 0.0;


for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){

 r = j*ngrids_X+i;

 if(i==0){
 first = psi[(j)*ngrids_X+i+1] + psi[(j)*ngrids_X + ngrids_X-1];
 }else if(i==(ngrids_X-1)){
 first = psi[(j)*ngrids_X] + psi[(j)*ngrids_X+i-1];
 }else{
 first = psi[(j)*ngrids_X+i+1] + psi[(j)*ngrids_X+i-1];
 }

 if(j==0){
 second = psi[(j+1)*ngrids_X+i] + psi[(ngrids_Y-1)*ngrids_X+i];
 }else if(j==(ngrids_Y-1)){
 second = psi[i] + psi[(j-1)*ngrids_X+i];
 }else{
 second = psi[(j+1)*ngrids_X+i] + psi[(j-1)*ngrids_X+i];
 }

 //psi1[r] = 0.50*(Sqr(dy)*first + Sqr(dx)*second + Sqr(dx*dy)*omega[r])*(1.0/(Sqr(dx) + Sqr(dy)));

 psi1[r] = (1 - sorcoeff)*psi[r] + sorcoeff*0.25*(first + second + Sqr(dx)*omega[r]);

 sumerr += Sqr(psi1[r] - psi[r]);

 }
}

for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){
 r = j*ngrids_X+i;
  
 psi[r] = psi1[r]; 

 }
}

 error = sqrt(sumerr/(ngrids_X*ngrids_Y));
 printf("%0.16lf\t%d\n",error,iter);
 iter++;
} 


for (i=0;i<ngrids_X;i++){
 for (j=0;j<ngrids_Y;j++){
 r = j*ngrids_X+i;
 
  if(j == 0){
  first = psi[ngrids_X + i] - psi[(ngrids_Y-1)*ngrids_X + i];
  //Uxp = v[(ngrids_Y-1)*ngrids_X + nx];
  }else if(j == ngrids_Y-1){
  first = psi[i] - psi[(ngrids_Y-2)*ngrids_X + i];
  //Uxp = v[(ngrids_Y-2)*ngrids_X + nx];
  }else{
  first = psi[(j+1)*ngrids_X + i] - psi[(j-1)*ngrids_X + i];
  //Uxp = v[(ny-1)*ngrids_X + nx];
  }

        if(i == 0){
	  second = psi[(j)*ngrids_X + 1] - psi[(j)*ngrids_X + ngrids_X-1];
	  //Uyp = u[(ny)*ngrids_X + ngrids_X-1];
	}else if(i == ngrids_X-1){
	  second = psi[(j)*ngrids_X] - psi[(j)*ngrids_X + ngrids_X-2];
	  //Uyp = u[(ny)*ngrids_X + ngrids_X-2];
	}else{
	  second = psi[(j)*ngrids_X + i+1] - psi[(j)*ngrids_X + i-1];
	  //Uyp = u[(ny)*ngrids_X + nx-1];
	}

  v[r] = (first/dy);
   
  u[r] = -(second/dx);

 }
}

vortex(v,u,dx,dy);


}//end of time


free(omega);
free(psi1);
//fftw_free(fomega);
free(psi);
//fftw_free(fpsi);
free(kx);
free(ky);
//free(fu);
//free(fv);
free(v);
free(u);
free(u1);
//free(v1);
free(in);
//free(out);
free(vorticity);
return 0;
}
