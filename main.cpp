//Librerias
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "loader.cpp"
#include "algorithms.c"
#include "main.h"

using namespace std;

/* Variables globales */
enum BOOL ILSA; /* Una variable booleana para saber si estamos ejecutando con el algoritmo de aceleración de bucle interno (ILA) 
o no. ILSA encuentra el bucle interno óptimo explorando todas las posibilidades. */
enum BOOL NOISOLATE;
enum BOOL USERDATA;
enum BOOL PARAMS;
enum BOOL LIMIT_DISTANCE;
#ifdef DYNALLOC
int LENGTH;
unsigned char *RNA1; 
unsigned char *RNA; /* Contains RNA string in terms of 0, 1, 2, 3 for A, C, G and U respectively*/
int *structure; /* An array to contain the optimal structure */
int *V; /* int V[LENGTH][LENGTH]; */
int *W;
int **VBI; /* VBI(i,j) will contain the energy of optimal internal loop closed with (i,j) base pair */
int **VM; /* VM(i, j) will contain the energy of optimla multiloop closed with (i,j) base pair */
int **WM; /* This array is introduced to help multiloop calculations. WM(i,j) contains the optimal energy of string segment from si to sj if this forms part of a multiloop */
int *indx; /* This array is used to index V array. Here V array is mapped from 2D to 1D and indx array is used to get the mapping back.*/
int *constraints;

double **QB;  // QB[i][j] is the sum over all possible loops closed by (i,j),
              // including the summed contributions of their subloops
double **Q;   // Q[i][j] in addition to the above quantity QB[i][j], Q[i][j]
              // also includes all configurations with (i,j) not paired
double **QM;  // QM[i][j] is the sum of configuration energies from i to j,
              // assuming that i,j are contained in a multiloop
double **P;   // P[i][j] The probability that nucleotides i and j form a basepair
#else
/* This are previously used variables, now they are not used. */
unsigned char RNA[LENGTH];
unsigned char RNA1[LENGTH];
int structure[LENGTH];
int VBI[LENGTH][LENGTH];
int VM[LENGTH][LENGTH];
int V[(LENGTH-1)*(LENGTH)/2 + 1]; /* int V[LENGTH][LENGTH]; */
int WM[LENGTH][LENGTH];
int W[LENGTH];
int indx [LENGTH];
#endif

// Funcion para calcular el tiempo
double segundos() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

// Inicialización de variables globales
void init_variables(int len) {
	int i;
#ifdef DYNALLOC
	LENGTH = len + 1;
	RNA = (unsigned char *) malloc(LENGTH * sizeof(unsigned char));
	RNA1 = (unsigned char *) malloc(LENGTH * sizeof(unsigned char));
	structure = (int *) malloc(LENGTH * sizeof(int));
	V = (int *) malloc(((LENGTH - 1) * (LENGTH) / 2 + 1) * sizeof(int));
	W = (int *) malloc(LENGTH * sizeof(int));
	VBI = (int **) malloc(LENGTH * sizeof(int *));
	for (i = 0; i < LENGTH; i++) {
		VBI[i] = (int *) malloc(LENGTH * sizeof(int));
	}
	VM = (int **) malloc(LENGTH * sizeof(int *));
	for (i = 0; i < LENGTH; i++) {
		VM[i] = (int *) malloc(LENGTH * sizeof(int));
	}
	WM = (int **) malloc(LENGTH * sizeof(int *));
	for (i = 0; i < LENGTH; i++) {
		WM[i] = (int *) malloc(LENGTH * sizeof(int));
	}
	indx = (int *) malloc(LENGTH * sizeof(int));
	constraints = (int*) malloc((len + 1) * sizeof(int));
#endif
	return;
}
/* Funcion principal
 *  1) Leer los argumentos de la línea de comandos.
 *  2) populate() de loader.cpp para leer los parámetros termodinámicos definidos 
 *     en los archivos dados en el directorio de datos.
 *  3) Inicializar variables
 *  4) Las llamadas que calculan la función definida en algoritmos.c para rellenar las tablas de energía.
 *  */
int main(int argc, char** argv) {
	int i;
	ifstream cf;
	int bases;
	string s, seq;
	int energy;
	float energia;
	double t1;
	ILSA = FALSE;
	NOISOLATE = FALSE;
	/* Leyendo linea de comandos */
	int fileIndex = 0, consIndex = 0, dataIndex = 0, paramsIndex=0, lcdIndex = 0;
	i = 1;
	while (i < argc) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i], "-ilsa") == 0) {
				ILSA = TRUE;
			} else if (strcmp(argv[i], "-noisolate") == 0) {
				NOISOLATE = TRUE;
			} else if (strcmp(argv[i], "-constraints") == 0) {
				if (i < argc)
					consIndex = ++i;
			} else if (strcmp(argv[i], "-params")==0) { 
				PARAMS = TRUE;			  
				if (i < argc)
					paramsIndex = ++i;
			} else if (strcmp(argv[i], "-datadir") == 0) {
				USERDATA = TRUE;
				if (i < argc)
					dataIndex = ++i;
			} else if (strcmp(argv[i], "-limitCD") == 0)
			{
				if (i < argc)
					lcdIndex = ++i;
			}
		} else {
			fileIndex = i;
		}
		i++;
	}
	cf.open(argv[fileIndex]);
	if (cf != NULL)
		printf("Archivo abierto.\n\n");
	else {
		printf("Error al abrir el archivo.\n\n");
		exit(-1);
	}
	seq = "";
	s = "";
	//Para archivo tipo FASTA
	char ss[10000];
	cf.getline(ss, 10000);
	if (ss[0] != '>') {
		char *fline;
		fline = strtok(ss, " ");
		while (fline != NULL) {
			seq.append(fline);
			fline = strtok(NULL, " ");
		}
	}
	while (!cf.eof()) {
		cf >> s;
		seq.append(s);
		s = "";
	}
	s = seq;
	bases = s.length();
	init_variables(bases);
	cout << "Secuencia ingresada: " << s << endl;   // Se imprime la secuencia
	cout << "Largo de la secuencia: " << bases <<endl; // Se imprime el largo de la secuencia
	cf.close();
	int **fbp = NULL, **pbp = NULL;
	if (handle_IUPAC_code(s, bases)  == FAILURE)
	{
		exit(0);
	}	
	if(USERDATA==TRUE)
		populate(argv[dataIndex],true);
	else if (PARAMS == TRUE)
		populate(argv[paramsIndex],false);
	else
		populate("combinaciones",false); //Lectura de archivos termodinámicos
	initTables(bases); // Se inicializan variables globales de acuerdo a las bases de la secuencia
	t1 = segundos();
	energy = calculate(bases, fbp, pbp, 0, 0); /* Ejecuta el algoritmo de programación dinámica para calcular
	 la energía óptima. Definido en el archivo algorithms.c*/
    t1 = segundos() - t1;
    energia = energy/100.00;
	cout << "\nEnergia minima libre: " << energia <<endl;
	cout << "\nEl calculo demoró "<< t1 << " segundos" <<endl;
	return 0;
}

GTFOLD_FLAGS handle_IUPAC_code(const string& s, const int bases)
{
	int* stack_unidentified_base;
	int stack_count=0;
	bool unspecd=0;
	stack_unidentified_base=new int[bases];
	// SH: Conversion de la secuencia en valores numericos.
	for(int i = 1; i <= bases; i++) {
		RNA[i] = getBase(s.substr(i-1,1));
		RNA1[i] = getBase1(s.substr(i-1,1));
		if (RNA[i]=='X') {
			return FAILURE; //exit(0);
		}
		else if(RNA[i]!='X' && RNA1[i]=='N'){
			unspecd=1;
			stack_unidentified_base[stack_count]=i;
			stack_count++;
		}
	}
}
