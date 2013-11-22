#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_matrix.h>


int main(int argc, char **argv){

FILE *input;
  int filas;
int numeros[10] = {0};
float *promi;
int columnas=24;
int s;
double *valor_propio;
  int i,j,k;
 int var;
  int grado;

  printf("This program is executed with %d arguments:\n", argc);
  for(i=0;i<argc;i++){
    printf("%s\n", argv[i]);
  }

  if(argc<2){
    printf("Se necesitan al menos 2 argumentos para ejecutar el programa!\n EXIT!\n");
    exit(1);
  }

 input=fopen(argv[1],"r");
    if(!input){
      printf("surgió un problema abriendo el archivo\n");
      exit(1);
    }

 filas=0;
    while ((var = fgetc(input)) != EOF){
      if (var =='\n')
	++filas;
    }





input=fopen(argv[1],"r");
gsl_matrix *DATOS = gsl_matrix_alloc (filas, columnas);
gsl_matrix_fscanf (input, DATOS);




 promi = malloc(filas * sizeof(float));
double a;
double prom;

for(i=0;i<filas;i++){
		 a=0.0;
		 prom=0.0;
	for(j=0;j<columnas;j++){
		double b=gsl_matrix_get (DATOS,i,j);
		a=a+b;
		}
	prom=a/columnas;
	promi[i]=prom;
	

	
}


gsl_matrix *COVA = gsl_matrix_alloc (filas, filas);
double suma;
for(i=0;i<filas;i++){

	for(j=0;j<filas;j++){
		suma=0.0;
		for(k=0;k<columnas;k++){

			double dik=gsl_matrix_get (DATOS,i,k);
			double djk=gsl_matrix_get (DATOS,j,k);
			double di=promi[i];
			double dj=promi[j];
			double e=(dik-di)*(djk-dj);
			suma=suma+e;
			}
	
	gsl_matrix_set (COVA, i, j, suma/columnas-1);

		}
}

	
valor_propio=malloc(filas * sizeof(double));
gsl_vector *eval = gsl_vector_alloc (filas);
gsl_matrix *evec = gsl_matrix_alloc (filas,filas);

gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (filas);
gsl_eigen_symmv (COVA, eval, evec, w);
  
gsl_eigen_symmv_free (w);
  
gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_ASC);

    for (i = 0; i < filas; i++)
      {
double eval_i= gsl_vector_get (eval, i);
gsl_vector_view evec_i= gsl_matrix_column (evec, i);
valor_propio[i]=eval_i;


gsl_vector_fprintf (stdout,&evec_i.vector, "%g");
      }


fclose(input);

input=fopen("OscarGarzon_eigenvectors.dat","w");

for(i=0;i<10;i++)
{
gsl_vector_view evec_i=gsl_matrix_column(evec,i);
gsl_vector_fprintf(input,&evec_i.vector,"%f ");
}
fclose(input);

input=fopen("OscarGarzon_eigenvalues.dat","w");
for(i=0;i<filas;i++)

{
fprintf(input, "%f\n", valor_propio[i]);

}







  }
  
  




