/* Generates a coordinates and velocities grid from trajectory data */
/* This grid is useful for plotting quiver plots */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"global.h"
#define	NDIM 2
#define Sqr(x) ((x)*(x))

/*

void vortex(double v[], double u[], double dy, double dx){

double Uxp, Uxn, Uyp, Uyn;
 //int nx, ny;
FILE *fpvort;
fpvort = fopen("vort_values","w");
//printf("hit void%f %lf\n",dx,dy);
//double Uxp, Uxn, Uyp, Uyn;

for(ny = 0 ; ny < ngrids_Y ; ny ++){
      for(nx = 0 ; nx < ngrids_X ; nx ++){
	n = (ny-1)*ngrids_X + nx;
	//printf("hi = %lf %lf",dx,dy);

}}



    for(ny = 1 ; ny <= ngrids_Y ; ny ++){
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny-1)*ngrids_X + nx;
	printf("hi = %lf %lf",v[n],u[n]);
	if(nx == 1){
	  Uyn = u[(ny-1)*ngrids_X + 2];
	  Uyp = u[(ny-1)*ngrids_X + ngrids_X];
	}else if(nx == ngrids_X){
	  Uyn = u[(ny-1)*ngrids_X + 1];
	  Uyp = u[(ny-1)*ngrids_X + ngrids_X-1];
	}else{
	  Uyn = u[(ny-1)*ngrids_X + nx+1];
	  Uyp = u[(ny-1)*ngrids_X + nx-1];
	}
	if(ny == 1){
	  Uxn = v[ngrids_X + nx];
	  Uxp = v[(ngrids_X-1)*ngrids_X + nx];
	   //Uxp = 0.0;
	}else if(ny == ngrids_Y){
	  Uxn = v[nx];
	  //Uxn = 0.0;
	  Uxp = v[(ngrids_Y-2)*ngrids_X + nx];
	}else{
	  Uxn = v[ny*ngrids_X + nx];
	  Uxp = v[(ny-2)*ngrids_X + nx];
	}
	
	vorticity[n] = (Uyn - Uyp)/(dx) - (Uxn - Uxp)/(dy);
	printf("vorticity = %lf", vorticity[n]);
      }
    }
fclose(fpvort);
}
*/






int main(int argc, char **argv){
  /*printf("%d\n",argc);
  if(argc != 5){
    printf("Usage : ./grid [*.xyz] [# of frames] [ngrids_X] [ngrids_Y]\n");
    exit(1);
  }*/

int iteration;

  int Frames = 400;
  char dummy[25];
  int frn;
  int nAtom, n, atomType, c;
  double timeNow, region[NDIM+1], regionH[NDIM+1]; 
  double *rx, *ry, *vx, *vy;
  double deltaT = 0.001;

  ngrids_X = 512;
  ngrids_Y = 512;
 
  int NHIST = 4;
  /*double **Grid;
  Grid = (double **)malloc((NHIST+1)*sizeof(double*));
  for(n = 0 ; n <= NHIST ; n ++)
    Grid[n] = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  int j;
  for(j = 1 ; j <= NHIST ; j ++)
    for(n = 1; n <= ngrids_Y*ngrids_X; n ++)
      Grid[j][n] = 0.;*/

  double *omega;
  omega = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  vorticity = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));


  
  /*double *streamFun;
  streamFun = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  
  double *Q;
  Q = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));*/

  double *vd;
  vd = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  double *ud;
  ud = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  //double *p;
  p = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  double *pd;
  pd = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *p1;
  p1 = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *density;
  density = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  //double *v;
  v = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  //double *u;
  u = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  double *velox;
  velox = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *veloy;
  veloy = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  
  double *rovx;
  rovx = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *rovy;
  rovy = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  double *ro;
  ro = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));

  rx = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));
  ry = (double *)malloc((ngrids_Y*ngrids_X+1)*sizeof(double));



  L = 10; 
  double dx, dy; 
  halfx = 0.5*L; 
  halfy = 0.5*L;
  dx = (double)L/(double)ngrids_X;
  dy = (double)L/(double)ngrids_Y;

  double deltax = 2.0*dx;
  double deltay = 2.0*dy;

  double gridWidth[NDIM+1];
 

//Reading of initial profile is done here


FILE *fpvort;
fpvort = fopen("QUIVER_time_0","r");

for(n = 1 ; n <= ngrids_X*ngrids_Y ; n ++){
   
  fscanf(fpvort,"%lf %lf %lf %lf %lf %lf", &rx[n], &ry[n], &vd[n], &ud[n], &omega[n], &rovx[n]);
  //printf("%lf %lf\n",vd[n], ud[n]);
  }
fclose(fpvort);





//printing of initial profile is done here
FILE *fpinitial1;
fpinitial1 = fopen("initial_0","w");

for(nx = 1 ; nx <= ngrids_X ; nx ++){
   
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){

    n = (ny - 1 )*ngrids_X + nx;
    x = (nx-0.5)*dx - halfx;
    y = (ny-0.5)*dy - halfy;
    fprintf(fpinitial1,"%lf %lf %lf %lf %lf\n", rx[n], ry[n], vd[n], ud[n], omega[n]);
  }
}

fclose(fpinitial1);


double mach = 0.2;
double U_0 = 0.25;
double Cs = U_0/mach;
double first, second, third, fourth, fifth, sixth;
double first1, second1, third1, fourth1, fifth1, sixth1;
double reynold = 0;
double viscosity = 0.01;


//Solving density equation

for(n = 1 ; n <= ngrids_X*ngrids_Y ; n ++){
ro[n] = 1.0;
}


int t;
//Start of time
  
for(t = 1 ; t <= 80000 ; t ++){





  for(nx = 1 ; nx <= ngrids_X ; nx ++){
   
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){

	n = (ny - 1 )*ngrids_X + nx;

      
	if(nx ==1){

	first = vd[(ny-1)*ngrids_X + 2] - vd[(ny-1)*ngrids_X + ngrids_X];
	second = ro[(ny-1)*ngrids_X + 2] - ro[(ny-1)*ngrids_X + ngrids_X];
	
	}else if(nx == ngrids_X){

	first = vd[(ny-1)*ngrids_X + 1] - vd[(ny-1)*ngrids_X + ngrids_X -1];
	second = ro[(ny-1)*ngrids_X + 1] - ro[(ny-1)*ngrids_X + ngrids_X-1];
	
	}else{

	first = vd[(ny-1)*ngrids_X + nx + 1] - vd[(ny-1)*ngrids_X + nx -1];
	second = ro[(ny-1)*ngrids_X + nx + 1] - ro[(ny-1)*ngrids_X + nx-1];
	}



	if(ny == 1){

	second1 = ro[ngrids_X + nx] - ro[(ngrids_Y-1)*ngrids_X + nx];
	first1 = ud[ngrids_X + nx] - ud[(ngrids_Y-1)*ngrids_X + nx];

	}else if(ny == ngrids_Y){

	second1 = ro[nx] - ro[(ngrids_Y-2)*ngrids_X + nx];
	first1 = ud[nx] - ud[(ngrids_Y-2)*ngrids_X + nx];

	}else{

	second1 = ro[(ny)*ngrids_X + nx] - ro[(ny-2)*ngrids_X + nx];
       	first1 = ud[(ny)*ngrids_X + nx] - ud[(ny-2)*ngrids_X + nx];

	}	

 
	ro[n] = ro[n]  - (ro[n]*first + ro[n]*first1 + vd[n]*second + ud[n]*second1);

     }
   }



  
  for(nx = 1 ; nx <= ngrids_X ; nx ++){
   
    for(ny = 1 ; ny <= ngrids_Y ; ny ++){

	n = (ny - 1 )*ngrids_X + nx;

	if(nx ==1){

	first = vd[(ny-1)*ngrids_X + 2] - vd[(ny-1)*ngrids_X + ngrids_X];
	fifth = ro[(ny-1)*ngrids_X + 2] - ro[(ny-1)*ngrids_X + ngrids_X];
	sixth = vd[(ny-1)*ngrids_X + 2] - 2.0*vd[(ny-1)*ngrids_X + 1] + vd[(ny-1)*ngrids_X + ngrids_X];
	first1 = ud[(ny-1)*ngrids_X + 2] - ud[(ny-1)*ngrids_X + ngrids_X];
	
	}else if(nx == ngrids_X){

	first = vd[(ny-1)*ngrids_X + 1] - vd[(ny-1)*ngrids_X + ngrids_X -1];
	fifth = ro[(ny-1)*ngrids_X + 1] - ro[(ny-1)*ngrids_X + ngrids_X-1];
	sixth = vd[(ny-1)*ngrids_X + 1] - 2.0*vd[(ny-1)*ngrids_X + ngrids_X] + vd[(ny-1)*ngrids_X + ngrids_X -1];
	first1 = ud[(ny-1)*ngrids_X + 1] - ud[(ny-1)*ngrids_X + ngrids_X -1];
	
	}else{

	first = vd[(ny-1)*ngrids_X + nx + 1] - vd[(ny-1)*ngrids_X + nx -1];
	fifth = ro[(ny-1)*ngrids_X + nx + 1] - ro[(ny-1)*ngrids_X + nx-1];
	sixth = vd[(ny-1)*ngrids_X + nx + 1] - 2.0*vd[(ny-1)*ngrids_X + nx] + vd[(ny-1)*ngrids_X + nx - 1];
	first1 = ud[(ny-1)*ngrids_X + nx + 1] - ud[(ny-1)*ngrids_X + nx -1];
	}

	if(ny == 1){

	second = vd[ngrids_X + nx] - vd[(ngrids_Y-1)*ngrids_X + nx];
	second1 = ud[ngrids_X + nx] - ud[(ngrids_Y-1)*ngrids_X + nx];
	fifth1 = ro[ngrids_X + nx] - ro[(ngrids_Y-1)*ngrids_X + nx];
	sixth1 = ud[ngrids_X + nx] - 2.0*ud[nx] + ud[(ngrids_Y-1)*ngrids_X + nx];	

	}else if(ny == ngrids_Y){

	second = vd[nx] - vd[(ngrids_Y-2)*ngrids_X + nx];
	second1 = ud[nx] - ud[(ngrids_Y-2)*ngrids_X + nx];
	fifth1 = ro[nx] - ro[(ngrids_Y-2)*ngrids_X + nx];
	sixth1 = ud[nx] - 2.0*ud[(ngrids_Y-1)*ngrids_X + nx] + ud[(ngrids_Y-2)*ngrids_X + nx];	

	}else{

	second = vd[(ny)*ngrids_X + nx] - vd[(ny-2)*ngrids_X + nx];
       	second1 = ud[(ny)*ngrids_X + nx] - ud[(ny-2)*ngrids_X + nx];
	fifth1 = ro[(ny)*ngrids_X + nx] - ro[(ny-2)*ngrids_X + nx];
	sixth1 = ud[(ny)*ngrids_X + nx] - 2.0*ud[(ny-1)*ngrids_X + nx] + ud[(ny-2)*ngrids_X + nx];	

	}	

        //vd[n] = vd[n] + (-vd[n]*(first/deltax) - ud[n]*(second/deltay)  - Sqr(Cs)*(fifth/(deltax*ro[n])) + (sixth/((Sqr(dx))*reynold*ro[n])))*deltaT;
        //ud[n] = ud[n] + (-vd[n]*(first1/deltax) - ud[n]*(second1/deltay) - Sqr(Cs)*(fifth1/(deltay*ro[n])) + (sixth1/((Sqr(dy))*reynold*ro[n])))*deltaT;
        vd[n] = vd[n] + (-vd[n]*(first/deltax) - ud[n]*(second/deltay)  - Sqr(Cs)*(fifth/(deltax*ro[n])) + viscosity*(sixth/((Sqr(dx))*ro[n])))*deltaT;
        ud[n] = ud[n] + (-vd[n]*(first1/deltax) - ud[n]*(second1/deltay) - Sqr(Cs)*(fifth1/(deltay*ro[n])) + viscosity*(sixth1/((Sqr(dy))*ro[n])))*deltaT;

 }
}




 
	

 double Uxp, Uxn, Uyp, Uyn;

    for(ny = 1 ; ny <= ngrids_Y ; ny ++){ 
      for(nx = 1 ; nx <= ngrids_X ; nx ++){
	n = (ny-1)*ngrids_X + nx;
	
	if(nx == 1){
	  Uyn = ud[(ny-1)*ngrids_X + 2];
	  Uyp = ud[(ny-1)*ngrids_X + ngrids_X];
	}else if(nx == ngrids_X){
	  Uyn = ud[(ny-1)*ngrids_X + 1];
	  Uyp = ud[(ny-1)*ngrids_X + ngrids_X-1];
	}else{
	  Uyn = ud[(ny-1)*ngrids_X + nx+1];
	  Uyp = ud[(ny-1)*ngrids_X + nx-1];
	}
	if(ny == 1){
	  Uxn = vd[ngrids_X + nx];
	  Uxp = vd[(ngrids_Y-1)*ngrids_X + nx];
	  //Uxp = 0.0;
	}else if(ny == ngrids_Y){
	  Uxn = vd[nx];
	  //Uxn = 0.0;
	  Uxp = vd[(ngrids_Y-2)*ngrids_X + nx];
	}else{
	  Uxn = vd[ny*ngrids_X + nx];
	  Uxp = vd[(ny-2)*ngrids_X + nx];
	}
	x = (nx-0.5)*dx - halfx;
	y = (ny-0.5)*dy - halfy;
	
	vorticity[n] = (Uyn - Uyp)/(2.0*dx) - (Uxn - Uxp)/(2.0*dy);
	//if((period == 1) || (T1 % 100) ==0){
	//if((T1 == (tfinal-1)*500))
	//if(t%1000 == 0)
	//fprintf(outfp1,"%lf\t %lf\t %lf\t %lf\t %lf\n", rx[n] , ry[n], vd[n], ud[n], vorticity[n]);
       //}
      }
    }

    if(t%100 == 0){

    FILE *outfp1;
    char outfile1[25];
    sprintf(outfile1,"Vorticity_states_%d", t);
    outfp1 = fopen(outfile1, "w");

    for(ny = 1 ; ny <= ngrids_Y ; ny ++){ 
      for(nx = 1 ; nx <= ngrids_X ; nx ++){

	n = (ny-1)*ngrids_X + nx;

	fprintf(outfp1,"%lf\t %lf\t %lf\t %lf\t %lf\n", rx[n] , ry[n], vd[n], ud[n], vorticity[n]);

      }
    } 	


    fclose(outfp1);

   }





}


 

  //free(Q);
  //free(omega);
  //free(streamFun);
  free(rx);
  free(ry);
  //free(vx);
  //free(vy);
  free(vd);
  free(p);
  free(p1);
  free(density);
  free(omega);
  free(v);
  free(u);
  free(velox);
  free(veloy);
  free(pd);
 
  return 0;
}
