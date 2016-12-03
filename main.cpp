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
// Algoritmos asociados
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
int LENGTH;
unsigned char *RNA1; 
unsigned char *RNA; // Contiene cadena de ARN en términos de 0, 1, 2, 3 para A, C, G y U respectivamente
int *V; // Para el calculo de MFE
int *W; // Para el calculo de MFE
int *constraints; // Ayuda al calculo de MFE
int **VBI; // VBI(i,j) contendrá la energía del bucle interno óptimo cerrado con el par de bases (i,j)
int **VM; // VM(i, j) contendrá la energía optima del multiloop cerrado con el par de bases (i,j)
int **WM; // Esta matriz se presenta para ayudar a los cálculos de multiloop. WM (i, j) contiene la energía óptima del segmento de cuerda de si a sj si esto forma parte de un multiloop
int *indx; // Esta matriz se utiliza para indexar V array. Aquí V matriz se asigna de 2D a 1D y indx matriz se utiliza para obtener la asignación de nuevo.
double **QB;  // QB [i] [j] es la suma de todos los bucles posibles cerrados por (i, j), incluyendo las contribuciones sumadas de sus subloops
double **Q;   // Q [i] [j] además de la cantidad QB [i] [j], Q [i] [j] incluye también todas las configuraciones con (i, j) no pareadas
double **QM;  // QM [i] [j] es la suma de las energías de configuración de i a j, suponiendo que i, j están contenidos en un multiloop

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
	//structure = (int *) malloc(LENGTH * sizeof(int));
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

BANDERA handle_IUPAC_code(const string& s, const int bases)
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
