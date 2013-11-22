#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 

int main(int argc, char **argv){
FILE *input;
  int num_lineas;
int s;
  int i,j;
 int var;
  int grado;
double chi_2;

  /*See what comes from the command line*/
  printf("This program is executed with %d arguments:\n", argc);
  for(i=0;i<argc;i++){
    printf("%s\n", argv[i]);
  }

  if(argc<3){
    printf("Se necesitan al menos 2 argumentos para ejecutar el programa!\n EXIT!\n");
    exit(1);
  }

 input=fopen(argv[1],"r");
    if(!input){
      printf("surgió un problema abriendo el archivo\n");
      exit(1);
    }

 grado=atoi(argv[2]);

printf("El grado del polinomio es %d \n", grado);

  if(grado<0)
    {
      printf("El NUmero que entra como parametro es Negativo! \nEXIT\n");
      exit(1);
    }
 


  num_lineas=0;
    while ((var = fgetc(input)) != EOF){
      if (var =='\n')
	++num_lineas;
    }

    printf("Número de lineas del archivo:\n -->%d\n",num_lineas); 

 float *x;
 float *y;
    x = malloc(num_lineas * sizeof(float));
    y = malloc(num_lineas * sizeof(float));
rewind(input);
    for(j=0;j<num_lineas;j++){
      fscanf(input, "%f %f",&x[j],&y[j]);

     } 

gsl_matrix *G = gsl_matrix_alloc (num_lineas, grado+1);
gsl_matrix *F = gsl_matrix_alloc (num_lineas,1);

 for (i = 0; i < num_lineas; i++){
    for (j = 0; j <= grado; j++){
 double h=pow(x[i],j);
      gsl_matrix_set (G, i, j, h);
}
}

for (i = 0; i < num_lineas; i++)
{

gsl_matrix_set (F,i,0,y[i]);


}

gsl_matrix *D= gsl_matrix_alloc (grado+1,grado+1);
 gsl_matrix_set_zero (D);
 gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                  1.0, G, G,
                  0.0, D);
gsl_permutation * p = gsl_permutation_alloc(grado+1);
gsl_matrix *invD = gsl_matrix_alloc(grado+1,grado+1);
gsl_linalg_LU_decomp (D, p, &s);
gsl_linalg_LU_invert(D,p,invD);

gsl_matrix *A= gsl_matrix_alloc (grado+1,num_lineas);
 gsl_matrix_set_zero (A);
gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0,invD, G, 0.0, A); 


gsl_matrix *RESPUESTA= gsl_matrix_alloc (grado+1,1);
gsl_matrix_set_zero (RESPUESTA);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,A, F, 0.0, RESPUESTA); 

printf("Los coeficientes del polinomio son \n");

for (i= 0; i<grado+1; i++)
{
double a=gsl_matrix_get (RESPUESTA,i,0);

printf("%f \n", a);}


double chi=0.0;
for(i=0;i<num_lineas;i++)
{
double fx=0.0;
for(j=0;j<grado+1;j++){

double a=gsl_matrix_get (RESPUESTA,j,0);
fx=fx+a*pow(x[i],j);

}

double b=gsl_matrix_get (F,i,0);
chi=pow(fx-b,2)+chi;



}

chi_2=chi/num_lineas;
printf("chi_2 =%f\n",chi_2);





}
/*---------------------------------------------------------
COmienzo de operaciones con matrices.
-----------------------------------------------------------*/


