//Librerias
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <math.h>
//#include "main-c.h"
// Auxiliares para el manejo de los archivos de datos
#define xstr(s) str(s)
#define str(s) #s
#define maxfil 100    /* maximum length of file names */
#define INFINITY_ 9999999  /* an arbitrary value given to infinity */
#define SMALLINFTY_ 99999 /*an arbitray value to determine safe range from infinity */
#define maxtloop 100 /* maximum tetraloops allowed (info read from tloop) */
#define maxstructures 1010 /* maximum number of structures in ct file */
#define maxbases 10000   /* maximum number of bases in a structure */
#define ctheaderlength 125 /* maximum length of string containing info on sequence */
#define ga_bonus -10 /* the value of a bonus for the "almost coaxial stacking" case in efn2 */
#define amax 400 /* this is a maximum line length for void linout (below) */
#define col 80  /* this is the number of columns in an output file */
#define numlen 8  /* maximum digits in a number */
#define maxforce 600 /* maximum number of bases that can be forced single */
#define maxgu 5 /* maximum number of u's in gu pair */
#define C_ 1 /* "c" for optimized VBI. */

#define MAXLOOP 30 /* Maximo tamaño de loop. */
#define MAXENG 1000

#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_U 3

using namespace std;

//Variables globales


// Esto era del main-c.h
#define DYNALLOC

#ifdef DYNALLOC
extern __thread int LENGTH;	
extern __thread unsigned char *RNA1; 
extern __thread unsigned char *RNA; // Contiene cadena de ARN en términos de 0, 1, 2, 3 para A, C, G y U respectivamente
extern __thread int *V; // Para el calculo de MFE
extern __thread int *W; // Para el calculo de MFE
extern __thread int *constraints; // Ayuda al calculo de MFE
extern __thread int **VBI; // VBI(i,j) contendrá la energía del bucle interno óptimo cerrado con el par de bases (i,j)
extern __thread int **VM; // VM(i, j) contendrá la energía optima del multiloop cerrado con el par de bases (i,j)
extern __thread int **WM; // Esta matriz se presenta para ayudar a los cálculos de multiloop. WM (i, j) contiene la energía óptima del segmento de cuerda de si a sj si esto forma parte de un multiloop
extern __thread int *indx; // Esta matriz se utiliza para indexar V array. Aquí V matriz se asigna de 2D a 1D y indx matriz se utiliza para obtener la asignación de nuevo.
#else

#define LENGTH 8500

extern unsigned char RNA[LENGTH];
extern unsigned char RNA1[LENGTH];
extern int structure[LENGTH];
extern int V[(LENGTH-1)*(LENGTH)/2 + 1];
extern int VBI[LENGTH][LENGTH];
extern int VM[LENGTH][LENGTH];
extern int WM[LENGTH][LENGTH];
extern int W[LENGTH];
extern int indx[LENGTH];
#endif

#define fourBaseIndex(a, b, c, d) (((a) << 6) + ((b) << 4) + ((c) << 2) + (d))
// Finaliza main-c.h

int poppen[5];
int maxpen;
int eparam[11]; // Ayuda al calculo de parametros similares
int multConst[3]; // Respaldo de las penalizaciones para multiloops
int dangle[4][4][4][2]; // Energias colindantes
int inter[31]; // Tamaño de penalización para loops internos
int bulge[31]; // Tamaño de penalización para bultos
int hairpin[31]; // Tamaño de penalizaciòn para loops de Horquill
int stack[256]; // Pilas de energía para pilas loop
int tstkh[256]; // Pilas de energia para loops de horquilla
int tstki[256]; //  Terminal de pilas de energia para loops internos
int tloop[maxtloop + 1][2];
int numoftloops; // Numero de loops
int iloop22[5][5][5][5][5][5][5][5]; // 2*1 loops internos
int iloop21[5][5][5][5][5][5][5]; // 2*1 loops internos
int iloop11[5][5][5][5][5][5]; // 1*1 loops internos
int coax[6][6][6][6]; 
int tstackcoax[6][6][6][6];
int coaxstack[6][6][6][6];
int tstack[6][6][6][6];
int tstkm[6][6][6][6];

int auend; // Para penalizacion AU
int gubonus;
int cint; // cint, cslope y c3 son usados por poly C loops de horquilla
int cslope;
int c3;
int efn2a; // Auxiliares del calculo de energia
int efn2b; // Auxiliares del calculo de energia
int efn2c; // Auxiliares del calculo de energia
int triloop[maxtloop + 1][2];
int numoftriloops;
int init;
int gail;
float prelog; // Usado por los loops que tengan tamaño mayor a 30

string EN_DATADIR;

//Funciones
void populate(const char *userdatadir,bool userdatalogic);
unsigned char getBase(string base);
unsigned char getBase1(string base);
int InicioValoresPila(string fileName);
int initMiscloopValues(string fileName);
int initDangleValues(string fileName);
int initLoopValues(string fileName);
int initTstkhValues(string fileName);
int initTstkiValues(string fileName);
int initTloopValues(string fileName);
int initInt21Values(string fileName);
int initInt22Values(string fileName);
int initInt11Values(string fileName);

//Para el manejo de los archivos de información
void populate(const char *userdatadir,bool userdatalogic) {

#ifndef GENBIN
	if (!userdatalogic) {
		EN_DATADIR.assign(xstr(DATADIR));
		EN_DATADIR += "/";
		EN_DATADIR += userdatadir;
	} else {
		EN_DATADIR.assign(userdatadir);
	}
#endif

	//Para el final del archivo
	if (EN_DATADIR[EN_DATADIR.length() - 1] != '/') {
		EN_DATADIR += "/";
	}
	initMiscloopValues("miscloop.dat"); // Bucles variados.
	InicioValoresPila("stack.dat"); // Energías libres para el apilamiento de pares de bases.
	initDangleValues("dangle.dat"); // Una sola base que apila energías libres.
	initLoopValues("loop.dat"); // Componente entropico para los bucles internos, del bulto y de horquilla.
	initTstkhValues("tstackh.dat"); // Free energies for terminal mismatch stacking in hairpin loops
	initTstkiValues("tstacki.dat"); // Energías libres para el empalme de la falta de coincidencia terminal en bucles horquilla
	initTloopValues("tloop.dat"); // Energías libres para tetraloops distinguidos
	initInt21Values("int21.dat"); // Energías libres para bucles interiores 2 x 1
	initInt22Values("int22.dat"); // Energías libres para bucles interiores 2 x 2
	initInt11Values("int11.dat"); // Energías libres para loops interiores 1 x 1 
}

// Transforma las bases a valores numericos
char baseADigito(string base) {
	if (!strcmp(base.c_str(), "A")) {
		return '1';
	}
	if (!strcmp(base.c_str(), "C")) {
		return '2';
	}
	if (!strcmp(base.c_str(), "G")) {
		return '3';
	}
	if (!strcmp(base.c_str(), "U")) {
		return '4';
	}
	if (!strcmp(base.c_str(), "N")) {
		return '5';
	}
	return (char) NULL;
}

unsigned char getBase(std::string base) {
	if (!strcmp(base.c_str(), "A") || !strcmp(base.c_str(), "a")) {
		return BASE_A;
	}
	if (!strcmp(base.c_str(), "C") || !strcmp(base.c_str(), "c")) {
		return BASE_C;
	}
	if (!strcmp(base.c_str(), "G") || !strcmp(base.c_str(), "g")) {
		return BASE_G;
	}
	if (!strcmp(base.c_str(), "U") || !strcmp(base.c_str(), "u") || !strcmp(
			base.c_str(), "T") || !strcmp(base.c_str(), "t")) {
		return BASE_U;
	}
	if (!strcmp(base.c_str(), "N") || !strcmp(base.c_str(), "n")||!strcmp(base.c_str(), "R")|| !strcmp(base.c_str(), "r")|| !strcmp(base.c_str(), "Y")|| !strcmp(base.c_str(), "y")|| !strcmp(base.c_str(), "M")|| !strcmp(base.c_str(), "m")|| !strcmp(base.c_str(), "K")|| !strcmp(base.c_str(), "k")|| !strcmp(base.c_str(), "S")|| !strcmp(base.c_str(), "s")|| !strcmp(base.c_str(), "W")|| !strcmp(base.c_str(), "w")|| !strcmp(base.c_str(), "B")|| !strcmp(base.c_str(), "b")|| !strcmp(base.c_str(), "D")|| !strcmp(base.c_str(), "d")|| !strcmp(base.c_str(), "H")|| !strcmp(base.c_str(), "h")|| !strcmp(base.c_str(), "V")|| !strcmp(base.c_str(), "v")) {
		return BASE_A;
	}
	return 'X';
}

// Para considerar distintas notaciones de cada base
unsigned char getBase1(std::string base) {
	if (!strcmp(base.c_str(), "A") || !strcmp(base.c_str(), "a")) {
		return BASE_A;
	}
	if (!strcmp(base.c_str(), "C") || !strcmp(base.c_str(), "c")) {
		return BASE_C;
	}
	if (!strcmp(base.c_str(), "G") || !strcmp(base.c_str(), "g")) {
		return BASE_G;
	}
	if (!strcmp(base.c_str(), "U") || !strcmp(base.c_str(), "u") || !strcmp(
			base.c_str(), "T") || !strcmp(base.c_str(), "t")) {
		return BASE_U;
	}
	if (!strcmp(base.c_str(), "N") || !strcmp(base.c_str(), "n")|| !strcmp(base.c_str(), "R")|| !strcmp(base.c_str(), "r")|| !strcmp(base.c_str(), "Y")|| !strcmp(base.c_str(), "y")|| !strcmp(base.c_str(), "M")|| !strcmp(base.c_str(), "m")|| !strcmp(base.c_str(), "K")|| !strcmp(base.c_str(), "k")|| !strcmp(base.c_str(), "S")|| !strcmp(base.c_str(), "s")|| !strcmp(base.c_str(), "W")|| !strcmp(base.c_str(), "w")|| !strcmp(base.c_str(), "B")|| !strcmp(base.c_str(), "b")|| !strcmp(base.c_str(), "D")|| !strcmp(base.c_str(), "d")|| !strcmp(base.c_str(), "H")|| !strcmp(base.c_str(), "h")|| !strcmp(base.c_str(), "V")|| !strcmp(base.c_str(), "v")) {
		return 'N';
	}
	return 'X';
}

int InicioValoresPila(string fileName) {

	ifstream cf; // Manejo de archivos
	int i, j, k, l; // Variables de apoyo
	int ii, jj, kk, ll; // Variables de apoyo al calculo de bases
	int index;
	char currentLine[256];
	string currentString;
	string s;
	// Inicializar el arreglo con infinito
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					stack[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "Error al abrir el archivo de consulta" << endl;
		exit(-1);
	}
	// Lectura de parametros termodinámicos, ayuda con el analisis de la secuencia en archivo
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}
	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;
	while (i < 16) {
		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);
		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;
		ll = 0;
		jj = 0;
		int z = 0;
		int r = 0;
		while (s[z] != '\0') {
			if (s[z] == ' ')
				z++;
			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			} else {
				char value[10];
				int x = 0;

				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}
				value[x] = '\0';
				int temp = (int) floor(100.0 * atof(value) + .5);
				stack[fourBaseIndex(ii,jj,kk,ll)] = temp;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;
		kk = (i % 4);
	}
	cf.close();
	return 0;
}

int initMiscloopValues(string fileName) {
	/*
	 miscloop.dat - Archivo de bucle. Contiene :
	 1. Extrapolación para bucles grandes basados ​​en la teoría de polímeros
	 2. Parámetros de corrección de bucle interno asimétricos
	 4. Paremetros para bucles multi-ramas (Solo para energia efn2)
	 6. Penalizaciones para combinaciones AU o GU
	 7. Bonus para horquilla GGG
	 */
	char currentWord[256];
	string s;
	ifstream cf; //Manejo de archivos
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "El archivo no se ha podido abrir" << endl;
		exit(-1);
	}
	s = "";
	cf >> currentWord;
	for (int index = 1; index < 13; index++) { // Lectura de 12 valores
		while (strcmp(currentWord, "-->")) {
			cf >> currentWord;
		}
		if (index == 1) {
			cf >> currentWord;
			prelog = 100 * atof(currentWord);
		}
		if (index == 2) {
			cf >> currentWord;
			maxpen = int(atof(currentWord) * 100.0 + .5);
		}
		if (index == 3) {
			for (int count = 1; count <= 4; count++) {
				cf >> currentWord;
				s = currentWord;
				poppen[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
		}
		if (index == 4) {
			eparam[1] = 0;
			eparam[2] = 0;
			eparam[3] = 0;
			eparam[4] = 0;
			eparam[7] = 30;
			eparam[8] = 30;
			eparam[9] = -500;
			int table[4];
			table[1] = 5;
			table[2] = 6;
			table[3] = 10;
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				multConst[count - 1] = (int) (atof(s.c_str()) * 100 + 0.5);
				eparam[table[count]] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
		}
		if (index == 5) {
			int table[4];
			for (int count = 1; count <= 3; count++) {
				cf >> currentWord;
				s = currentWord;
				table[count] = (int) (atof(s.c_str()) * 100 + 0.5);
			}
			efn2a = table[1];
			efn2b = table[2] - 1;
			efn2c = table[3] - 1;
		}
		if (index == 6) {
			cf >> currentWord;
			auend = (int) (100 * atof(currentWord));
		}
		if (index == 7) {
			cf >> currentWord;
			gubonus = (int) (100 * atof(currentWord));
		}
		if (index == 8) {
			cf >> currentWord;
			cslope = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 9) {
			cf >> currentWord;
			cint = (int) (100 * atof(currentWord));
		}
		if (index == 10) {
			cf >> currentWord;
			c3 = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 11) {
			cf >> currentWord;
			init = (int) (100 * atof(currentWord)) + 1;
		}
		if (index == 12) {
			cf >> currentWord;
			gail = (int) floor(.5 + atof(currentWord));
		}
	}
	cf.close();
	return 0;
}

int initDangleValues(string fileName) {
	ifstream cf; // Manejo de archivos
	char currentLine[256];
	string currentString;
	string s;
	int index;
	int i, j, k, l;
	int ii, jj, kk, ll; // ii = 1st base, jj = 2nd base, kk = 3rd base, ll = 1 arriba o 2 abajo
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 2; l++) {
					dangle[i][j][k][l] = INFINITY_;
				}
			}
		}
	}
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "La apertura del archivo a fallado" << endl;
		exit(-1);
	}

	// The 8 first lines are junk
	for (index = 1; index <= 8; index++) {
		cf.getline(currentLine, 256);
	}

	// 8 lines of useful data
	i = 0;
	ii = 0;
	jj = 0;
	kk = 0;
	ll = 0;

	while (i < 8) {

		if (i != 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);

		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;

		jj = 0;

		int z = 0;
		int r = 0;

		while (s[z] != '\0') {

			if (s[z] == ' ')
				z++;

			else if (s[z] == '.') {
				z++;
				kk++;
				if (kk == 4)
					kk = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			} else {
				char value[10];
				int x = 0;
				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}
				value[x] = '\0';
				int temp = (int) floor(100.0 * atof(value) + .5);
				dangle[ii][jj][kk][ll] = temp;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				kk++;
				if (kk == 4)
					kk = 0;
			}
		}
		i++;
		ii++;
		if (ii == 4)
			ii = 0;
		if (i == 4)
			ll = 1;
	}
	cf.close();
	return 0;
}

int initLoopValues(string fileName) {
	// algorithm.c, linea 2996
	ifstream cf; // Archivo actual
	char currentLine[256];
	char currentWord[256];
	string s;
	int index;
	int tempValue = 0;
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "Error al abrir el archivo" << endl;
		exit(-1);
	}
	// Para la basura del archivo
	for (index = 1; index <= 4; index++) {
		cf.getline(currentLine, 256);
	}
	while (index < 30) {
		for (int j = 1; j <= 4; j++) {
			cf >> currentWord;
			if (j == 1)
				index = atoi(currentWord);
			if (j > 1) {
				if (strcmp(currentWord, ".")) {
					tempValue = (int) (100 * atof(currentWord) + 0.5);
				} else {
					tempValue = INFINITY_;
				}
			}
			switch (j) {
			case 2:
				inter[index] = tempValue;
				break;
			case 3:
				bulge[index] = tempValue;
				break;
			case 4:
				hairpin[index] = tempValue;
				break;
			}
		}
	}
	cf.close();
	return 0;
}

int initTstkhValues(string fileName) {
	ifstream cf; //cf = Archivo actual
	int i, j, k, l;
	int ii, jj, kk, ll;
	int index;
	char currentLine[256];
	string currentString;
	string s;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					tstkh[fourBaseIndex(i,j,k,l)] = INFINITY_; // Calculo para la matriz
				}
			}
		}
	}
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "Error al abrir el archivo" << endl;
		exit(-1);
	}
	// Para basura del archivo
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}
	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;
	while (i < 16) {
		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);
		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;
		ll = 0;
		jj = 0;
		int z = 0;
		int r = 0;
		while (s[z] != '\0') {
			if (s[z] == ' ')
				z++;
			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			}

			else {
				char value[10];
				int x = 0;
				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}
				value[x] = '\0';

				int temp = (int) floor(100.0 * atof(value) + .5);
				tstkh[fourBaseIndex(ii,jj,kk,ll)] = temp;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;
		jj = 0;
		kk = (i % 4);
	}
	cf.close();
	return 0;
}

int initTstkiValues(string fileName) {
	ifstream cf; //cf = current file
	int i, j, k, l;
	int ii, jj, kk, ll;
	int index;
	char currentLine[256];
	string currentString;
	string s;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				for (l = 0; l < 4; l++) {
					tstki[fourBaseIndex(i,j,k,l)] = INFINITY_;
				}
			}
		}
	}
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "Error al abrir el archivo" << endl;
		exit(-1);
	}
	for (index = 1; index <= 15; index++) {
		cf.getline(currentLine, 256);
	}

	i = 0;
	kk = 0;
	ii = 0;
	jj = 0;
	ll = 0;
	while (i < 16) {
		if (i % 4 == 0)
			for (index = 1; index < 9; index++)
				cf.getline(currentLine, 256);
		cf.getline(currentLine, 256);
		s = currentLine;
		j = 0;
		ll = 0;
		jj = 0;
		int z = 0;
		int r = 0;
		while (s[z] != '\0') {
			if (s[z] == ' ')
				z++;
			else if (s[z] == '.') {
				z++;
				ll++;
				if (ll == 4)
					ll = 0;
				r++;
				if (r % 4 == 0)
					jj++;
			}
			else {
				char value[10];
				int x = 0;
				while (s[z] != ' ' && s[z] != '\0') {
					value[x++] = s[z++];
				}
				value[x] = '\0';
				int temp = (int) floor(100.0 * atof(value) + .5);
				tstki[fourBaseIndex(ii,jj,kk,ll)] = temp;
				r++;
				z++;
				if (r % 4 == 0)
					jj++;
				ll++;
				if (ll == 4)
					ll = 0;
			}
		}
		i++;
		if (!(i % 4))
			ii++;
		jj = 0;
		kk = (i % 4);
	}
	cf.close();
	return 0;
}

// Inicializa los valores para los tipos de loops
int initTloopValues(string fileName) {
	ifstream cf;
	int count;
	char currentLine[256];
	char currentSeqNumbers[7];
	char currentValue[6];

	currentSeqNumbers[6] = '\0';
	currentValue[5] = '\0';

	string s, temp;
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		cerr << "Error al abrir el archivo" << endl;
		exit(-1);
	}
	cf.getline(currentLine, 256);
	cf.getline(currentLine, 256);
	numoftloops = 0;
	while (!cf.eof() && (++(numoftloops) < maxtloop)) {
		int clindex=0;
		cf.getline(currentLine, 256);
		while(currentLine[clindex]== ' ') clindex++;
		for (count = 0; count < 6; count++) {
			temp = currentLine[count + clindex];
			currentSeqNumbers[count] = baseADigito(temp);
		}
		clindex=clindex+7;
		while(currentLine[clindex]== ' ') clindex++;
		count = 0;
		while(currentLine[clindex+count]!=' '&&currentLine[clindex+count]!='\0') {
			currentValue[count] = currentLine[count + clindex];
			count++;
		}
		tloop[numoftloops][0] = (int) atoi(currentSeqNumbers);
		tloop[numoftloops][1] = (int) floor(100.0 * atof(currentValue) + 0.5);
	}
	cf.close();
	return 0;
}

int initInt22Values(string fileName) {

	int i, j, k, r, q, t, y, z;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							for (y = 0; y < 4; y++)
								for (z = 0; z < 4; z++)
									iloop22[i][j][k][r][q][t][y][z] = INFINITY_;
	ifstream cf;
	int index, flag;
	char currentLine[256], currentValue[6];
	string s, s1, s2, temp;
	string sre;
	int base[4];
	int l, m;
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) {
		exit(-1);
	}
	for (index = 1; index < 28; index++) {
		cf.getline(currentLine, 256);
	}
	sre = "Y";
	flag = 0;
	for (index = 1; index <= 36; index++) { // 36 tablas 16x16
		// Principio de las tablas
		while (!flag) {
			cf.getline(currentLine, 256);
			s = currentLine;
			int z = 0;
			while (s[z] != '\0') {
				if (s[z] == 'Y')
					flag = 1;
				z++;
			}
		}
		flag = 0;
		// Saltamos 5 lineas
		for (i = 0; i < 5; ++i) {
			cf.getline(currentLine, 256);
		}
		// Obtener las bases de cierre
		cf.getline(currentLine, 256);
		s1 = currentLine;
		cf.getline(currentLine, 256);
		s2 = currentLine;
		s = s1 + s2;
		int z = 0;
		int k = 0;
		while (s[z] != '\0') {
			if (s[z] == 'A')
				base[k++] = BASE_A;
			else if (s[z] == 'C')
				base[k++] = BASE_C;
			else if (s[z] == 'G')
				base[k++] = BASE_G;
			else if (s[z] == 'U')
				base[k++] = BASE_U;
			z++;
		}
		cf.getline(currentLine, 256);
		for (int rowIndex = 1; rowIndex <= 16; rowIndex++) {
			for (int colIndex = 1; colIndex <= 16; colIndex++) {
				cf >> currentValue;
				j = ((rowIndex - 1) - (rowIndex - 1) % 4) / 4;
				k = (rowIndex - 1) % 4;
				l = ((colIndex - 1) - (colIndex - 1) % 4) / 4;
				m = (colIndex - 1) % 4;
				iloop22[base[0]][base[1]][base[2]][base[3]][j][l][k][m] = (int) floor(100.0 * atof(currentValue) + 0.5);
			}
		}
	}
	cf.close();
	return 0;
}

int initInt21Values(string fileName) {

	ifstream cf;
	char currentLine[256];
	string sre;
	string s, s1, s2;
	int a, b, c, d, e, f, g, index;
	int i, j, k, r, q, t, y;
	int z;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							for (y = 0; y < 4; y++)
								iloop21[i][j][k][r][q][t][y] = INFINITY_;
	a = 0;
	b = 0;
	c = 0;
	d = 0;
	e = 0;
	f = 0;
	g = 0;
	k = 0;
	int base1[7];
	int base2[7];
	base1[1] = BASE_A + 1;
	base2[1] = BASE_U + 1;
	base1[2] = BASE_C + 1;
	base2[2] = BASE_G + 1;
	base1[3] = BASE_G + 1;
	base2[3] = BASE_C + 1;
	base1[4] = BASE_U + 1;
	base2[4] = BASE_A + 1;
	base1[5] = BASE_G + 1;
	base2[5] = BASE_U + 1;
	base1[6] = BASE_U + 1;
	base2[6] = BASE_G + 1;
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail())
		exit(-1);
	for (index = 1; index <= 17; index++) {
		cf.getline(currentLine, 256);
	}
	i = 1;
	while (i <= 6) {
		j = 1;
		while (j <= 4) {
			k = 1;
			for (index = 1; index <= 10; index++)
				cf.getline(currentLine, 256);
			s = currentLine;
			while (k <= 4) {
				cf.getline(currentLine, 256);
				s = currentLine;
				r = 0;
				z = 0;
				int jj = 1;
				d = 1;
				while (s[z] != '\0') {
					if (s[z] == ' ')
						z++;
					else if (s[z] == '.') {
						z++;
						d++;
						if (d == 5)
							d = 1;
						r++;
						if (r % 4 == 0)
							jj++;
					}
					else {
						char value[10];
						int x = 0;
						while (s[z] != ' ' && s[z] != '\0') {
							value[x++] = s[z++];
						}
						value[x] = '\0';
						int temp = (int) floor(100.0 * atof(value) + .5);
						a = base1[i];
						b = base2[i];
						f = base1[jj];
						g = base2[jj];
						c = k;
						e = j;
						iloop21[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1][g - 1] = temp;
						r++;
						z++;
						if (r % 4 == 0)
							jj++;
						d++;
						if (d == 5)
							d = 1;
					}
				}
				k++;
			}
			j++;
		}
		i++;
	}
	cf.close();
	return 0;
}

int initInt11Values(string fileName) {

	int i, j, k, r, q, t;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				for (r = 0; r < 4; r++)
					for (q = 0; q < 4; q++)
						for (t = 0; t < 4; t++)
							iloop11[i][j][k][r][q][t] = INFINITY_;
	ifstream cf;
	int index;
	char currentLine[256];
	string s;
	int base1[7];
	int base2[7];
	int a, b, c, d, f;
	fileName = EN_DATADIR + fileName;
	cf.open(fileName.c_str(), ios::in);
	if (cf.fail()) 
		exit(-1);
	for (index = 1; index <= 17; index++) {
		cf.getline(currentLine, 256);
	}
	// Arreglos de 6x6: AU CG GC UA GU UG
	base1[1] = BASE_A + 1;
	base2[1] = BASE_U + 1;
	base1[2] = BASE_C + 1;
	base2[2] = BASE_G + 1;
	base1[3] = BASE_G + 1;
	base2[3] = BASE_C + 1;
	base1[4] = BASE_U + 1;
	base2[4] = BASE_A + 1;
	base1[5] = BASE_G + 1;
	base2[5] = BASE_U + 1;
	base1[6] = BASE_U + 1;
	base2[6] = BASE_G + 1;
	i = 0;
	k = 0;
	while (k < 6) {
		k++;
		index = 0;
		for (index = 1; index <= 10; index++) {
			cf.getline(currentLine, 256);
		}
		i = 0;
		b = 1;
		while (i < 4) {
			int jj = 1;
			++i;
			cf.getline(currentLine, 256);
			s = currentLine;
			j = 0;
			int r = 0;
			int z = 0;
			int e = 1;
			while (s[z] != '\0') {
				if (s[z] == ' ')
					z++;
				else if (s[z] == '.') {
					z++;
					e++;
					if (e == 5)
						e = 1;
					r++;
					if (r % 4 == 0)
						jj++;
				}
				else {
					char value[10];
					int x = 0;
					while (s[z] != ' ' && s[z] != '\0') {
						value[x++] = s[z++];
					}
					value[x] = '\0';
					int temp = (int) floor(100.0 * atof(value) + .5);
					a = base1[k];
					d = base2[k];
					c = base1[jj];
					f = base2[jj];
					iloop11[a - 1][b - 1][c - 1][d - 1][e - 1][f - 1] = temp;
					r++;
					z++;
					if (r % 4 == 0)
						jj++;
					e++;
					if (e == 5)
						e = 1;
				}
			}
			b++;
		}
	}
	cf.close();
	return 0;
}
