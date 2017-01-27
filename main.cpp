//Librerias
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Algoritmos asociados para el calculo
#include "loader.cpp"
#include "algorithms.c"
//#include "main.h"
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
__thread int LENGTH;	
__thread unsigned char *RNA1; 
__thread unsigned char *RNA; // Contiene cadena de ARN en términos de 0, 1, 2, 3 para A, C, G y U respectivamente
__thread int *V; // Para el calculo de MFE
__thread int *W; // Para el calculo de MFE
__thread int *constraints; // Ayuda al calculo de MFE
__thread int **VBI; // VBI(i,j) contendrá la energía del bucle interno óptimo cerrado con el par de bases (i,j)
__thread int **VM; // VM(i, j) contendrá la energía optima del multiloop cerrado con el par de bases (i,j)
__thread int **WM; // Esta matriz se presenta para ayudar a los cálculos de multiloop. WM (i, j) contiene la energía óptima del segmento de cuerda de si a sj si esto forma parte de un multiloop
__thread int *indx; // Esta matriz se utiliza para indexar V array. Aquí V matriz se asigna de 2D a 1D y indx matriz se utiliza para obtener la asignación de nuevo./*
__thread ifstream **cf;
/*
+ El problema principal en la paralelizacion con pthread en este algoritmo, es que estas variables globales que se 
+ encuentran arriba son modificadas en todas las funciones que son parte de la funcion de calculo del valor de energía
+ esto quiere decir que al trabajar con threads, todos accederán y modificarán las mismas variables lo que lleva a error
+ Se penso 
+
*/

// Funciones
enum BANDERA { SUCCESS = 0, FAILURE, ERR_OPEN_FILE, NO_CONS_FOUND};
BANDERA initialize_constraints(int*** fbp, int*** pbp, int& numpConstraints, int& numfConstraints, const char* constr_file);
BANDERA handle_IUPAC_code(const string& s, const int bases);
void limit_contact_distance(int lCD, int length);
void force_noncanonical_basepair(const char* nc_basepairs, int length);

string s;
int largoseq=0;
int numthreads;
//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; 

// Struct para el uso de N Posix Threads
struct Message {

   int start, stop, cadena, thread;
   
};

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

// Generacion de las secuencias aleatorias

double sre_random(void)
{
  static long  rnd1;    /* Numero random de LCG1 */
  static long  rnd2;    /* Numero random de LCG2 */
  static long  rnd;     /* Numero random de return */
  static long  tbl[64]; /* Tabla para Bays/Durham shuffle */
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
      for (i = 0; i < 64; i++) {
		  x = a1*(rnd1%q1);
		  y = r1*(rnd1/q1);
		  rnd1 = x-y;
		  if (rnd1 < 0)
			rnd1 += m1;
		  x = a2*(rnd2%q2); 
		  y = r2*(rnd2/q2);
		  rnd2 = x-y;
		  if (rnd2 < 0)
			rnd2 += m2;
		  tbl[i] = rnd1-rnd2;
		  if (tbl[i] < 0)
			tbl[i] += m1;
      }
      sre_randseed = 0;
  }
  x    = a1*(rnd1%q1);
  y    = r1*(rnd1/q1);
  rnd1 = x-y;
  if (rnd1 < 0)
	rnd1 += m1;
  x    = a2*(rnd2%q2);
  y    = r2*(rnd2/q2);
  rnd2 = x-y;
  if (rnd2 < 0) 
	rnd2 += m2;
  i   = (int) (((double) rnd / (double) m1) * 64.);
  rnd = tbl[i];
  tbl[i] = rnd1-rnd2;
  if (tbl[i] < 0)
	tbl[i] += m1;
  return ((double) rnd / (double) m1);  
}

// La proxima generación de calculos aleatorios
//Existe una serie de funciones de la libreria squid que realizan este proceso, esta funcion strshuffle es la mas simple
//pero todas tienen la misma funcion, la diferencia es el metodo que utilizan (una de ellas utiliza cadenas de markov)
int StrShuffle(string &s1, string s2)
{
  int  len;
  int  pos;
  char c;
  
  if (s1 != s2) s1=s2;
  for (len = s1.length(); len > 1; len--)
    {				
      pos       = CHOOSE(len);
      c         = s1[pos];
      s1[pos]   = s1[len-1];
      s1[len-1] = c;

    }
  return 1;
}


/* Funcion Paralela*/

void *Funcion(void *ptr){
	int p;
	string seq, s_random, segmento, copiaseg;
	int energy;
	float energia;
	double energy_random_prom,valor_extraido;
	double energy_random_desv[1000];
	double desv=0;
	double Z;
	//double t1;
	struct Message *data;
    data = (struct Message *) ptr;
	for(p = data->start; p < data->stop; p = p + 1){
		segmento = s.substr(data->cadena, largoseq);
		if((int)segmento.length() < largoseq){
			segmento = s.substr(data->cadena);
		}
		//cout <<"\nPara "<< p << " comienzo "<<data->start << " y final "<<data->stop<<" cadena: "<<data->cadena<<endl;	
		//pthread_mutex_lock(&mutex);
		init_variables(largoseq);
		int **fbp = NULL, **pbp = NULL;
		int numfConstraints = 0, numpConstraints = 0;
		if (handle_IUPAC_code(segmento, largoseq)  == FAILURE)
			exit(0);
		//populate("combinaciones",false); //Lectura de archivos termodinámicos
		initTables(largoseq); // Se inicializan variables globales de acuerdo a las bases de la secuencia
		energy = calculate(largoseq, fbp, pbp, 0, 0); 
	    //pthread_mutex_unlock(&mutex);
	    energia = energy/100.00;
	    //cout << "\nEnergia minima libre: " << energia<<" para segmento "<< segmento <<endl;
		for (int i=0; i<1000; i++){
			StrShuffle(copiaseg,segmento); // NO demora, no es necesario paralelizar
		//	pthread_mutex_lock(&mutex);
			init_variables(largoseq);
			if (handle_IUPAC_code(copiaseg, largoseq)  == FAILURE) // Para el error
				exit(0);
		//	populate("combinaciones",false);			
			initTables(largoseq);
			valor_extraido = calculate(largoseq, fbp, pbp, numfConstraints, numpConstraints);
		//	pthread_mutex_unlock(&mutex);
			energy_random_prom = valor_extraido + energy_random_prom; 
			energy_random_desv[i] = valor_extraido/100.00;
		}
		energy_random_prom = energy_random_prom/(100000);
		for (int j=0; j<1000; j++)
			desv = (energy_random_desv[j] - energy_random_prom)*(energy_random_desv[j] - energy_random_prom)+desv; 
		desv = sqrt(desv/1000);
		Z = (energia - (energy_random_prom))/desv;
		//printf("Thread: %i  Valor-Z = %f\n",data->thread, Z);
		cout << "En el segmento " << p  <<" : "<< segmento << " se tiene el valor Z " << Z <<endl;
		//cout <<"\nPara "<< p << " comienzo "<<data->start << " y final "<<data->stop<<" cadena: "<<data->cadena<<endl;
		data->cadena = data->cadena + largoseq/2;
		desv=0;
		Z=0;																																															
		energy_random_prom=0;
	}
	//}
	//t1 = segundos() - t1; // Calculo de tiempo para la energia de la secuencia completa
	//cout << "El calculo del Z-Score para todos los segmentos demoro "<<t1<<" segundos"<<endl<<endl;
	pthread_exit(0);
}

/* Funcion principal
 *  1) Leer los argumentos de la línea de comandos.
 *  2) populate() de loader.cpp para leer los parámetros termodinámicos definidos 
 *     en los archivos dados en el directorio de datos.
 *  3) Inicializar variables
 *  4) Las llamadas que calculan la función definida en algoritmos.c para rellenar las tablas de energía.
 *  */
int main(int argc, char** argv) {
	int i,h,fileIndex = 0;;
	double t1;
	ifstream cf;
	string seq, s_random, segmento, copiaseg;
	//double t1;
	cout <<  "Ingrese el largo de los segmentos: \n" ; 
	cin >> largoseq;
	cout <<  "Ingrese numero de threads: \n" ; 
	cin>> numthreads;

	ILSA = FALSE;
	NOISOLATE = FALSE;
	/* Leyendo linea de comandos */
	i = 1;
	while (i < argc) {
		if (argv[i][0] == '-') {
			/*if (i < argc){
				lcdIndex = ++i;			
			}*/
		} 
		else {
			fileIndex = i;
		}
		i++;
	}
	cf.open(argv[fileIndex]);
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
	//cout << "Secuencia ingresada: " << s << endl;   // Se imprime la secuencia
	//cout << "Largo de la secuencia: " << bases <<endl; // Se imprime el largo de la secuencia
	cf.close();
	//t1 = 0;
	//t1 = segundos();  // Se empieza a calcular el tiempo
/* Paralelo*/
	float NUM, redondeo;
	int NUM2, NUM3;//, NUM4;
	int secuencia=0, corredor=0,contador=0;
	secuencia = s.length();
	pthread_attr_t attribute;
	pthread_attr_init(&attribute);
    pthread_attr_setdetachstate(&attribute,PTHREAD_CREATE_JOINABLE);
    pthread_t *thread;
    void *exit_status;  
    thread = static_cast<pthread_t*>((void *)calloc(numthreads,sizeof(pthread_t)));
    struct Message **data;
    data = static_cast<struct Message**>(calloc(numthreads,sizeof(struct Message *))); 
    for (i = 0; i < numthreads; i = i + 1)
        data[i] = static_cast<struct Message*>(calloc(1,sizeof(struct Message)));
	while (corredor+largoseq<=secuencia){
		contador++;
		corredor = corredor + largoseq/2;
	}
	cout << "La cantidad de segmentos a analizar es: "<<contador<<endl;
	//NUM4 = secuencia/numthreads;
	NUM3 = contador;
	NUM2 = contador/numthreads;
	NUM = (float)(contador)/(float)(numthreads);
	redondeo = NUM-(float)(NUM2);
	if(redondeo>0.49)
		NUM2++;
	//int cadena = (largoseq/2) * NUM2;
	t1 = 0;
	t1 = segundos();  // Se empieza a calcular el tiempo para el Z-Score
	for(h = 0; h < numthreads; h = h + 1){
		if(h == 0){
			data[h]->start = 0;
			data[h]->stop = NUM2; 
			data[h]->cadena = 0;
			data[h]->thread = h;
		} 
		else{
			data[h]->thread = h;
			data[h]->cadena = ((largoseq/2) * NUM2) + data[h-1]->cadena;
			if(h == (numthreads - 1)){
				data[h]->start = data[h-1]->stop;
				data[h]->stop = NUM3;
//				populate("combinaciones",false);
			}
			else{
				data[h]->start = data[h-1]->stop;
				data[h]->stop = NUM2 + data[h-1]->stop;
//				populate("combinaciones",false);				
			}
		}
		populate("combinaciones",false);
		pthread_create(&thread[h], &attribute, Funcion, (void *) data[h]);
	}
	pthread_attr_destroy(&attribute);
	t1 = segundos() - t1; // Calculo de tiempo para la energia de la secuencia completa
	cout << "El calculo del Z-Score para todos los segmentos demoro "<<t1<<" segundos"<<endl<<endl;
    for (i = 0; i < numthreads; i = i + 1)
        pthread_join(thread[i],&exit_status);
	return 0;
}

BANDERA handle_IUPAC_code(const string& s, const int bases)
{
	int* stack_unidentified_base;
	int stack_count=0;
	//bool unspecd=0;
	stack_unidentified_base=new int[bases];
	// SH: Conversion de la secuencia en valores numericos.
	for(int i = 1; i <= bases; i++) {
		RNA[i] = getBase(s.substr(i-1,1));
		RNA1[i] = getBase1(s.substr(i-1,1));
		if (RNA[i]=='X') {
			return FAILURE; //exit(0);
		}
		else if(RNA[i]!='X' && RNA1[i]=='N'){
	//		unspecd=1;
			stack_unidentified_base[stack_count]=i;
			stack_count++;
		}
	}
}

bool is_valid_base(char c)
{	
	return ( (c-'A' == 0) || (c-'a' == 0) || 
			 (c-'C' == 0) || (c-'c' == 0) ||
			 (c-'G' == 0) || (c-'g' == 0) ||
			 (c-'U' == 0) || (c-'u' == 0));
}

void trim_spaces(string& str)
{
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );

}

void tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}
