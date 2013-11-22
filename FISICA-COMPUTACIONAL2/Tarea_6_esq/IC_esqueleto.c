/*

Este programa toma cuatro numeros como prametros.
Todos son del tipo float.
Genera un archivo de texto con las condiciones iniciales necesaria para modelar una galaxia simplificada.

*/
	

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
	
/* Se definen las variables y las constantes*/
	
	FILE *input;
	double G;
	double M;
	double PARSEC;
	double m;
	double x,y,v_x,v_y;
	int i,j;


int main(int argc, char **argv){

	
	/* Verificamos que el numero de parametros es el correcto */
	x=atof(argv[1]);
	y=atof(argv[2]);
	v_x=atof(argv[3]);
	v_y=atof(argv[4]);
	printf("Este programa se esta ejecutando con %d argumento(s):\n",argc-1);
	for(i=1;i<argc;i++){
    	printf("%s\n", argv[i]);
  	}
  	if(argc-1!=4){
   	printf("Se necesitan 4 parametros para determinar las condiciones iniciales\n EXIT!\n");
    	exit(1);
    	}

	/* Aqui irá un loop que utilizara la funcion asignarCondiciones(int orbita, double angulo) para dar a cada particula 
	su posicion y velocidad inicial*/

	/* Finalmente creamos un archivo y escribimos en él las condiciones de cada particula en una tabla*/

	input=fopen("Condiciones_iniciales.dat","w");
	
}

double asignarCondiciones(int orbita, double angulo){

	/*Método que asigna a una particula su velocidad y posicion inicial.
	Parametros: 
	int orbita: Determina la orbita a la que esta la particula en cuention. Con este parametro se determinan tambien la velocidad y7 ep 		radio.
	int angulo: Determina la posicion en la orbita de la particula.
	return: Velocidad y posicion de la particula medida con respecto a la masa central de la galaxia.
	*/
return 0;
}

