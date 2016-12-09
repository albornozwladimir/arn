//Librerias
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Algoritmos asociados para el calculo
#include "loader.cpp"
#include "algorithms.c"
#include "main.h"
// Para la generación de secuencias random
static int sre_randseed = 42;
#define CHOOSE(a)   ((int) (sre_random() * (a)))
//
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

/* deallocate global variables */
void free_variables() {
	int i;

#ifdef DYNALLOC
	free(indx);
	for (i = 0; i < LENGTH; i++)
		free(WM[i]);
	free(WM);
	for (i = 0; i < LENGTH; i++)
		free(VM[i]);
	free(VM);
	for (i = 0; i < LENGTH; i++)
		free(VBI[i]);
	free(VBI);
	free(W);
	free(V);
	free(constraints);
	free(RNA);
	free(RNA1);

#endif

	return;

}

double
sre_random(void)
{
  static long  rnd1;    /* random number from LCG1 */
  static long  rnd2;            /* random number from LCG2 */
  static long  rnd;             /* random number we return */
  static long  tbl[64];   /* table for Bays/Durham shuffle */
  long x,y;
  int i;

  /* Magic numbers a1,m1, a2,m2 from L'Ecuyer, for 2 LCGs.
   * q,r derive from them (q=m/a, r=m%a) and are needed for Schrage's algorithm.
   */
  long a1 = 40014;    
  long m1 = 2147483563;   
  long q1 = 53668;
  long r1 = 12211;

  long a2 = 40692;
  long m2 = 2147483399;
  long q2 = 52774;
  long r2 = 3791;

  if (sre_randseed > 0) 
    {
      rnd1 = sre_randseed;
      rnd2 = sre_randseed;
        /* Fill the table for Bays/Durham */
      for (i = 0; i < 64; i++) {
  x    = a1*(rnd1%q1);   /* LCG1 in action... */
  y    = r1*(rnd1/q1);
  rnd1 = x-y;
  if (rnd1 < 0) rnd1 += m1;

  x    = a2*(rnd2%q2);   /* LCG2 in action... */
  y    = r2*(rnd2/q2);
  rnd2 = x-y;
  if (rnd2 < 0) rnd2 += m2;

  tbl[i] = rnd1-rnd2;
  if (tbl[i] < 0) tbl[i] += m1;
      }
      sre_randseed = 0;   /* drop the flag. */
    }/* end of initialization*/


  x    = a1*(rnd1%q1);   /* LCG1 in action... */
  y    = r1*(rnd1/q1);
  rnd1 = x-y;
  if (rnd1 < 0) rnd1 += m1;

  x    = a2*(rnd2%q2);   /* LCG2 in action... */
  y    = r2*(rnd2/q2);
  rnd2 = x-y;
  if (rnd2 < 0) rnd2 += m2;

        /* Choose our random number from the table... */
  i   = (int) (((double) rnd / (double) m1) * 64.);
  rnd = tbl[i];
      /* and replace with a new number by L'Ecuyer. */
  tbl[i] = rnd1-rnd2;
  if (tbl[i] < 0) tbl[i] += m1;

  return ((double) rnd / (double) m1);  
}

int StrShuffle(string &s1, string s2)
{
  int  len;
  int  pos;
  int i=0;
  char c;
  
  if (s1 != s2) s1=s2;
  for (len = s1.length(); len > 1; len--)
    {				
      pos       = CHOOSE(len);
      c         = s1[pos];
      s1[pos]   = s1[len-1];
      s1[len-1] = c;

    }
    for(i; i <= s2.length(); i++){
        printf("%c", s1[i]);
      }
    printf("\n");
  return 1;
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
	string s, seq, s_random;
	int energy;
	float energia;
	double energy_random_prom;
	double energy_random_desv[50];
	double desv=0;
	double Z;
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
	int numfConstraints = 0, numpConstraints = 0;
	if (consIndex != 0)
	{
		GTFOLD_FLAGS r = initialize_constraints(&fbp, &pbp, numpConstraints, numfConstraints, argv[consIndex]);
		if (r == ERR_OPEN_FILE)
		{
			free_variables();
			exit(-1);
		}
	}
	
	
	if (handle_IUPAC_code(s, bases)  == FAILURE)
	{
		free_variables();
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
	energy = calculate(bases, fbp, pbp, numfConstraints, numpConstraints); /* Ejecuta el algoritmo de programación dinámica para calcular
	 la energía óptima. Definido en el archivo algorithms.c*/
    t1 = segundos() - t1;
    energia = energy/100.00;
	cout << "\nEnergia minima libre: " << energia <<endl;

	for(int i=0;i<5;i++){
    	StrShuffle(s,seq);
    	int bases2 = s.length();
    	init_variables(bases2);
    	int **fbp = NULL, **pbp = NULL;
		int numfConstraints = 0, numpConstraints = 0;		
		if (handle_IUPAC_code(s, bases2)  == FAILURE)
		{
			free_variables();
			exit(0);
		}		
		if(USERDATA==TRUE)
			populate(argv[dataIndex],true);
		else if (PARAMS == TRUE)
			populate(argv[paramsIndex],false);
		else
			populate("Turner99",false); /* Defined in loader.cc file to read in the thermodynamic parameter values from the tables in the ../data directory. */
    	initTables(bases2);
    	energy_random_prom = (calculate(bases2, fbp, pbp, numfConstraints, numpConstraints) + energy_random_prom);
    	printf("Valores-prom = %f\n", (calculate(bases2, fbp, pbp, numfConstraints, numpConstraints))/100.00);
    	energy_random_desv[i] =calculate(bases2, fbp, pbp, numfConstraints, numpConstraints)/100.00;
    
  	}
  	energy_random_prom = energy_random_prom/(5*100);
  	for(i=0;i<5;i++){
  		desv= (energy_random_desv[i] - energy_random_prom)*(energy_random_desv[i] - energy_random_prom)+desv;
  	}
  	printf("Valor-prom = %f\n", energy_random_prom);
  	desv = sqrt(desv/4);
  	printf("Valor-desv = %f\n", desv);
  	Z = (energia - (energy_random_prom))/desv;
  	printf("Valor-Z = %f\n", Z);
	cout << "\nEl calculo demoró "<< t1 << " segundos" <<endl;
	return 0;
}

GTFOLD_FLAGS initialize_constraints(int*** fbp, int ***pbp, int& numpConstraints, int& numfConstraints, const char* constr_file)
{
	ifstream cfcons;

	fprintf(stdout, "Running with constraints\n");
	//fprintf(stdout, "Opening constraint file: %s\n", argv[consIndex]);
	fprintf(stdout, "Opening constraint file: %s\n", constr_file);

	cfcons.open(constr_file, ios::in);
	if (cfcons != NULL)
		fprintf(stdout, "Constraint file opened.\n");
	else {
		fprintf(stderr, "Error opening constraint file\n\n");
		cfcons.close();
		return ERR_OPEN_FILE; //exit(-1);
	}

	char cons[100];

	while (!cfcons.eof()) {
		cfcons.getline(cons, 100);
		if (cons[0] == 'F' || cons[0] == 'f')
			numfConstraints++;
		if (cons[0] == 'P' || cons[0] == 'p')
			numpConstraints++;
	}
	cfcons.close();

	fprintf(stdout, "Number of Constraints given: %d\n\n", numfConstraints
			+ numpConstraints);
	if (numfConstraints + numpConstraints != 0)
		fprintf(stdout, "Reading Constraints.\n");
	else {
		fprintf(stderr, "No Constraints found.\n\n");
		return NO_CONS_FOUND;
	}

	*fbp = (int**) malloc(numfConstraints * sizeof(int*));
	*pbp = (int**) malloc(numpConstraints * sizeof(int*));

	int fit = 0, pit = 0, it = 0;

	for (it = 0; it < numfConstraints; it++) {
		(*fbp)[it] = (int*) malloc(2* sizeof (int));
	}
	for(it=0; it<numpConstraints; it++) {
		(*pbp)[it] = (int*)malloc(2*sizeof(int));
	}
	cfcons.open(constr_file, ios::in);

	while(!cfcons.eof()) {
		cfcons.getline(cons,100);
		char *p=strtok(cons, " ");
		p = strtok(NULL, " ");
		if(cons[0]=='F' || cons[0]=='f') {
			int fit1=0;
			while(p!=NULL) {
				(*fbp)[fit][fit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			fit++;
		}
		if( cons[0]=='P' || cons[0]=='p') {
			int pit1=0;
			while(p!=NULL) {
				(*pbp)[pit][pit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			pit++;
		}
	}

	fprintf(stdout, "Forced base pairs: ");
	for(it=0; it<numfConstraints; it++) {
		for(int k=1;k<= (*fbp)[it][2];k++)
			fprintf(stdout, "(%d,%d) ", (*fbp)[it][0]+k-1, (*fbp)[it][1]-k+1);
	}
	fprintf(stdout, "\nProhibited base pairs: ");
	for(it=0; it<numpConstraints; it++) {
		for(int k=1;k<= (*pbp)[it][2];k++)
			fprintf(stdout, "(%d,%d) ", (*pbp)[it][0]+k-1, (*pbp)[it][1]-k+1);
	}
	fprintf(stdout, "\n\n");
	
	return SUCCESS;
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
