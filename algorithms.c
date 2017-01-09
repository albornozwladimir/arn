#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data.h"
#include "constants.h"
#include "main-c.h"
#include "algorithms.h"
#include <omp.h> // Permite trabajar con computacion paralela
#define DEBUG 0
#define WM(i,j) WM[j][i]  /* Se definen los tipos que contendrán los calculos.*/

unsigned int chPairKey;

//Constantes
int plen = 0, flen = 0, sslen = 0;
int *pbpi, *pbpj, *fbpi, *fbpj, *ss;

// Esta función calcula que chPairKey sea procesado por la función chPair. 
void init_chPair() {
	int i, j;

	chPairKey = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			chPairKey += checkPair(i, j) << ((i << 2) + j);
}

int update_chPair(int i, int j) 
{
	int r = 0;
	if (!((i >= 0 && i <=3 )&&(j >=0 && j <=3)))
		return r;

	if (!(chPairKey & (1 << ((i << 2) + j))))
	{
		chPairKey += 1 << ((i << 2) + j);	
		r = 1;
	}
	return r;
}


// Este pragma devuelve 1 si la base b1 y b2 pueden emparejarse, de lo contrario devuelve 0, usando chPairKey calculado en la función init_chPair. Aquí b1 yb2 son 0-3 para representar uno de los cuatro nucleótidos A, C, G y U
#if 0
inline
int chPair(int b1, int b2) {
	return (chPairKey & (1 << ((b1<<2) + b2)));
}
#else
#define chPair(a, b)  (chPairKey & (1 << (((a)<<2) + (b))))  /* Please try to run this, to understand this statement. Defined by Professor Bader. */
#endif

// Inicialización de variables
void initTables(int len) {

	int i, j;
	int LLL;

#if 0
	int z = (len)*(len+1)/2 + 1;

	V = new int[z];
	indx = new int[len+1];
#endif

#if DEBUG
#ifdef DYNALLOC
	if (len != LENGTH-1)
		printf("ERROR: en initTables, largo (%5d) != Largo-1 (%5d)\n",len,LENGTH-1);
#endif
#endif

	init_chPair();

	for (i = 0; i < LENGTH; i++) {
		W[i] = INFINITY_; // Inicializando la matriz W con INFINITY
		constraints[i] = 0;
#if 0
		indx[i] = (LENGTH-1)*(i-1) - (i*(i-1))/2;
		indx[i] = (len)*(i-1) - (i*(i-1))/2;
#endif
		for (j = 0; j < LENGTH; j++) {
			VBI[i][j] = INFINITY_;
			VM[i][j] = INFINITY_;
			WM[i][j] = INFINITY_;
		}
	}
	LLL = (LENGTH - 1) * (LENGTH) / 2 + 1;
	for (i = 0; i < LLL; i++)
		V[i] = INFINITY_;
	for (i = 0; i <= LENGTH - 1; i++)
		indx[i] = (len) * (i - 1) - (i * (i - 1)) / 2;
	return;
}

int checkSS(int i, int j) {
	int it;
	for (it = i + 1; it < j; it++) {
		if (constraints[it] > 0)
			return 1;
	}
	return 0;
}

//Funcion principal que permite ejecutar los calculos de energia
int calculate(int len, int **forceList, int **prohibitList, int forcelen, int prohibitlen) {
	int b, i, j, it, k;
	for(i=1;i<=len;i++) 
	{
		if(RNA1[i]=='N') 
			constraints[i] = -1;
	}	
	if (prohibitlen != 0) 
	{
			for (it = 0; it < prohibitlen; it++) 
			{
				for(k= 1; k <= prohibitList[it][2];k++)
				{
					constraints[prohibitList[it][0]+k-1] = -1;
					if(prohibitList[it][1]!=0)
					{
						constraints[prohibitList[it][1]+1-k] = -1;
					}
				}
			}
	}
	if (forcelen != 0) 
	{
			for (it = 0; it < forcelen; it++) 
			{
				for(k=1; k <= forceList[it][2];k++)
				{
					if (!chPair(RNA[forceList[it][0]+k-1], RNA[forceList[it][1]-k+1])) 
					{
						continue;
					}
					constraints[forceList[it][0]+k-1] = forceList[it][1]+1-k;
					constraints[forceList[it][1]+1-k] = forceList[it][0]+k-1;
				}
			}
	}
	// Aquí b-1 es la longitud del segmento cerrado con (i, j) par de bases. Suponemos que el tamaño mínimo de un lazo horquilla cerrado con (i, j) igual a 3
	// Para b = 4 a 6, los bucles horquilla y en b = 6 bucles de la pila son posibles. Por lo tanto, sólo WM, y V matriz deben ser calculados.
	 // Si (i, j) no puede emparejarse entonces sólo necesita ser calculado WM.
	for (b = 4; b <= 6; b++) {
	#pragma omp parallel for private (i,j) schedule(guided)
		for (i = 1; i <= len - b; i++) {
			j = i + b;
			if (chPair(RNA[i], RNA[j])) // Comprueba si las bases i y j se juntan o no
				calcVWM(i, j, INFINITY_, INFINITY_); // Calcula el array V y WM para el elemento (i, j)
			else
				calcWM(i, j); // Calcula el array WM para el elemento (i, j)
		}
	}
	// Tener en cuenta que los cálculos de bucles internos utilizando el algoritmo de aceleración tiene que hacerse para cada par de bases de cierre (i, j) incluso si no es capaz de emparejarse.
    // * Para ocuparse de esto, ambos casos han sido separados usando una variable booleana ILSA
	if (ILSA == FALSE) { /* If we are executing internal loop speedup algorithm (ILSA) */
		/* For b=7 to 10, base pair (i,j) is not able to form multiloops. */
		for (b = 7; b <= 10; b++) {
		#pragma omp parallel for private (i,j) schedule(guided)
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				if (chPair(RNA[i], RNA[j])) {
					calcVBI(i, j); /* Calculates VBI element at (i,j) */
					calcVWM(i, j, VBI[i][j], INFINITY_); /* Calculates V and WM arrays*/
				} else
					calcWM(i, j); /* Calculates WM element at (i,j) */
			}
		}
		for (b = 11; b <= len - 1; b++) {
		#pragma omp parallel for private (i,j) schedule(guided)			
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				if (chPair(RNA[i], RNA[j])) {
					calcVBIVMVWM(i, j); // Calcula los elementos VBI, VM, V y WM en (i, j)
				} else
					calcWM(i, j); // Calcula el elemento WM en (i, j)
			}
		}
	} else { // Si estamos ejecutando con ILSA - algoritmo de aceleración de bucle interno
		for (b = 7; b <= 10; b++) {
		#pragma omp parallel for private (i,j) schedule(guided)
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				calcVBIS(i, j); // Calcula la matriz VBI [i] [j] con el algoritmo de aceleración de bucle interno (ILSA)
				if (chPair(RNA[i], RNA[j])) {
					calcVWM(i, j, VBI[i][j], INFINITY_); // Calcula el elemento V y WM en (i, j)
				} else {
					calcWM(i, j);
				}
			}
		}
		for (b = 11; b <= len - 1; b++) {
		#pragma omp parallel for private (i,j) schedule(guided)
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				calcVBIS(i, j); // Cálculo de la matriz VBI en (i, j) - Hecho en ambos casos si (i, j) se agrupa o no
				if (chPair(RNA[i], RNA[j])) {
					calcVMVWM(i, j); // Cálculo de WM, V, WM en orden en (i, j)
				} else
					calcWM(i, j); // Cálculo de WM en (i, j)
			}
		}
	}         
		for (j = 5; j <= len; j++) // La relación de recursión para la matriz W no depende de ninguna otra matriz, por lo que se puede hacer después de que el cálculo de otras matrices estén terminadas.
			calcW(j);
	return W[len];
}   // Fin de la funcion que calcula

/* Esta función calcula la energía óptima de bucles internos cerrados con pares de bases (i, j) usando una heurística, lo que limita su tamaño a un valor constante - MAXLOOP
 Un bucle interno contiene un par de bases de cierre (i, j) y un par de bases cerrado (ip, jp). Esta función busca el mejor par de bases cerrado para el par de bases de cierre (i, j) dentro de la ventana dada limitada por MAXLOOP
 */
void calcVBI(int i, int j) {
	int ip, jp, temp, VBIij, thres1;
	VBIij = INFINITY_;
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		return;
	// Teniendo ip = i + 1 y jp = j-1, crea un bucle de pila. Los bucles de pila se cuidan por separado en el cálculo de V utilizando la función eS (). Por lo tanto, para ip = i + 1, el valor jp debe ser menor o igual que j-2
	ip = i + 1;
	thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); 
	/*Minimo tamaño del loop de horquilla que puede cerrar el par (ip, jp) es 3 que resulta en el valor minimo de jp=ip+4  */
	for (jp = thres1; jp <= j - 2; jp++) {
		if (chPair(RNA[ip], RNA[jp])) {
			if (checkSS(i, ip) || checkSS(jp, j))
				continue;
			temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; 
			/* Energia del loop interno cerrado por (i,j) y (ip, jp) + la energía optima de la subestructura cerrada por (ip,jp)*/
			if (VBIij > temp)
				VBIij = temp;
		}
	}
	for (ip = i + 2; ip <= i + MAXLOOP + 1; ip++) {
		thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); 
		/* El tamaño minimo de un loop de horquilla es 3, entonces jp inicia de ip+4 */
		for (jp = thres1; jp <= j - 1; jp++) {
			if (chPair(RNA[ip], RNA[jp])) {
				if (checkSS(i, ip) || checkSS(jp, j))
					continue;
				temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; 
				/*Energía de un loop interno cerrado por (i,j) y (ip,jp) + la energía optima de la subestructura cerrada por (ip,jp)*/				
				if (VBIij > temp)
					VBIij = temp;
			}
		}
	}
	VBI[i][j] = VBIij;
	return;
}

/* Algoritmo de speedup para loops internos */
/* Calculo de loops internos utilizando este algoritmo de speedup. El algoritmo calcula el loop optimo cerrado por el par base (i,j) */
void calcVBIS(int i, int j) {
	int ip, jp, E, VBIij, c = 3, b, len = LENGTH - 1, E1, E2, g; /* ip and jp form enclosed base pairs and c is a small constant currently taken as 3. The loops having one or both sides smaller than c are calculated as special cases. */
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		return;
	VBIij = VBI[i][j];
	/*Caso 1: los loops tienen el primer lado mas corto que c y el segundo lado tiene todos los tamaños permitidos*/
	/*Teniendo ip=i+1 y jp=j-1 crea un stack el cual se trabaja con la funcion eS. Entonces aqui el valor maximo de jp será j-2 */
	ip = i + 1;
	for (jp = ip + 4; jp <= j - 2; jp++) {
		if (chPair(RNA[ip], RNA[jp])) {
			E = eL(i, j, ip, jp) + V[indx[ip] + jp];
			if (VBIij > E)
				VBIij = E;
		}
	}
	for (ip = i + 2; ip <= i + c; ip++) {
		for (jp = ip + 4; jp <= j - 1; jp++) { /* tamaño minimo de un loop de horquilla es 3.*/
			if (chPair(RNA[ip], RNA[jp])) {
				E = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > E)
					VBIij = E;
			}
		}
	}
	/*Caso 2: Cuando el primer lado es mayor o igual a c pero el segundo es mas pequeño*/
	for (ip = i + c + 1; ip < j - 1; ip++) {
		for (jp = j - c; jp <= j - 1 && jp >= ip + 4; jp++) { /* tamaño minimo de un loop de horquilla es 3.*/
			if (chPair(RNA[ip], RNA[jp])) {
				E = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > E)
					VBIij = E;
			}
		}
	}
	/*Caso 3: Caso general, cuando ambos lados de los loops internos son mayores o iguales a c*/
	/*Casos base para este (i,j) son g=j-i-2c-3 y j-i-2c-4, los valores de brecha deben ser siempre mayores o iguales a 3 */
	/* Primer caso base, ambos lados del loop son iguales a c*/
	ip = i + c + 1;
	jp = j - c - 1;
	g = jp - ip - 1;
	if (g < 3) {
		VBI[i][j] = VBIij;
		return;
	} /* si g es menor que 3, entonces no se extiende este valor. En este caso el segundo caso base g-1=j-1-2c-4 tampoco hará una valor de brecha válido*/
	E = eL(i, j, ip, jp) + V[indx[ip] + jp];
	if (VBIij > E)
		VBIij = E;
	/* Se extiende este caso base para todos los pares base de la forma (i-b, j+b)*/
	for (b = 1; b <= MIN(i - 1, len - j); b++) {
		E = eL(i - b, j + b, ip, jp) + V[indx[ip] + jp];
		/* Dos opciones mas para el par base (i-b,j+b), estas son introducidas teniendo uno de los lados el resultante loop interno exactamente igual a c */
		/* El segundo lado es c con el par base (i-b,j+b) cerrando  */
		int ip1 = i + c + 1 + b;
		int jp1 = (j + b) - c - 1;
		E1 = eL(i - b, j + b, ip1, jp1) + V[indx[ip1] + jp1];
		/* El primer lado es c con el par base (i-b,j+b) cerrando */
		int ip2 = (i - b) + c + 1;
		int jp2 = j - c - 1 - b;
		E2 = eL(i - b, j + b, ip2, jp2) + V[indx[ip2] + jp2];
		if (E > E1) {
			E = E1;
			ip = ip1;
			jp = jp1;
		}
		if (E > E2) {
			E = E2;
			ip = ip2;
			jp = jp2;
		}
		if (VBI[i - b][j + b] > E) {
			VBI[i - b][j + b] = E;
		}
	}
	if (g == 3) {
		VBI[i][j] = VBIij;
		return;
	} /* En este caso la brecha g-1=2, la cual no debiera ser extendida*/
	ip = i + c + 2;
	jp = j - c - 1;
	E1 = eL(i, j, ip, jp) + V[indx[ip] + jp];
	E2 = eL(i, j, i + c + 1, j - c - 2) + V[indx[i + c + 1] + j - c - 2];
	if (VBIij > E1)
		VBIij = E1;
	if (VBIij > E2)
		VBIij = E2;
	if (E2 < E1) {
		ip = i + c + 1;
		jp = j - c - 2;
	}
	for (b = 1; b <= MIN(i - 1, len - j); b++) {
		E = eL(i - b, j + b, ip, jp) + V[indx[ip] + jp];
		/*Dos mas opciones para el par base (i-b,j+b), teniendo uno de los lados iguales a c  */
		/* El primer lado es igual a c */
		int ip1 = (i - b) + c + 1;
		int jp1 = (j) - c - 2 - b;
		E1 = eL(i - b, j + b, ip1, jp1) + V[indx[ip1] + jp1];
		/* El segundo lado es igual a c */
		int ip2 = i + c + 2 + b;
		int jp2 = (j + b) - c - 1;
		E2 = eL(i - b, j + b, ip2, jp2) + V[indx[ip2] + jp2];
		if (E > E1) {
			E = E1;
			ip = ip1;
			jp = jp1;
		}
		if (E > E2) {
			E = E2;
			ip = ip2;
			jp = jp2;
		}
		if (VBI[i - b][j + b] > E) {
			VBI[i - b][j + b] = E;
		}
	}
	VBI[i][j] = VBIij;
}

/* Función para calcular el valor de WM(i,j) */
void calcWM(int i, int j) {
	int b = multConst[2], c = multConst[1]; /* b es la penalidad de brecha y c es la penalidad para bases singulares de multiloops */
	int h;
	/* WMidjd= base colgada en los lados i-ésimo y j-ésimo. 
	WMidj=base colgada en el lado i-ésimo.
	WMijd=base colgada en el lado j-ésimo.
	WMij = sin base colgada en ambos lados
	 */
	int WMidjd, WMidj, WMijd, WMij, WMijp;
	int rnai, rnaj;
	rnai = RNA[i]; /* Se lee el valor de RNA[i] y RNA[j] para hacer más rápida la ejecución del programa */
	rnaj = RNA[j];
	WMijp = INFINITY_;
	/* El tamaño minimo de un loop de horquilla es 3, eso hace que el limite de inicio sea h=i+4 y el de final sea j-5  */ 
	for (h = i + 4; h < j - 4; h++) {
		int temp = WM[i][h] + WM(h+1,j);
		if (temp <= WMijp)
			WMijp = temp;
	}
	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;
	/* Si la base i y j se parean*/
	WMij = V[indx[i] + j] + auPen(rnai, rnaj) + b;
	/* Si la base i+1 y j se parean. Se agrega la energia de la interacción colgante del par base (i+1,j) con la base i estando en el final 3'  */
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(RNA[i + 1], rnaj) + b + c;
	/* Si la base i y j-1 se parean. Se agrega la energia de la interacción colgante del par base (i,j-1) con la base j estando en el final 5'  */
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(rnai, RNA[j - 1]) + b + c;
	/* Si la base i+1 y j-1 se parean. Se agrega la energia de la interacción colgante del par base (i+1,j-1) con la base i estando en el final 3' y la base j estando en el final 5'  */
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1] + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1], RNA[j - 1]) + b + 2* c ;
	/* Se toma el minimo de todos los términos */
	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));
	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;
	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];
	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];
	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WMij = MIN(WMijp, WMij);
	WM[i][j] = WMij;
	WM(i,j) = WMij; 
	/* con esta instrucción se hace el arreglo WM simétrico. Esta macro convertirá esta instrucción en WM[j][i] = WM[i][j] haciendo que WM sea simétrico.  */
	return;
}

/* Función usada para calcular el valor de V y WM para un par i,j dado. El calculo de WM en (i,j) requiere el valor de V para (i,j) */
void calcVWM(int i, int j, int VBIij, int VMij) {
	int a, b, c, h, Vij, eh, es;
	int WMidjd, WMidj, WMijd, WMij, WMijp;
	int rnai, rnaj;
	rnai = RNA[i];
	rnaj = RNA[j];
	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;
	/* Se inicia V */
	eh = eH(i, j); /*Energía de un loop de horquilla */
	es = eS(i, j); /*Energía de un stack, con (i,j) y (i-1,j+1) pares de bases */
	if (es == 0) {
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];
	Vij = MIN(MIN(eh, es), MIN(VBIij, VMij));
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		Vij = INFINITY_;
	V[indx[i] + j] = Vij;
	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		//Se verifica si i+1,j-1 se han pareado
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			
			//Si no se revisa i-1,j+1
			//Con pares de bases isolados se verifica adelante
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			if (Vijl > INFINITY_ - SMALLINFTY_)
				//Si se encuentra un par de bases, la energía se setea a infinito
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		printf("ERROR: En calcVMW: i: %5d  j: %5d\n",i,j);
#endif
	/* V finaliza */
	/* WM inicia */
	a = multConst[0];
	b = multConst[2];
	c = multConst[1];
	WMijp = INFINITY_;
	for (h = i + 4; h <= j - 5; h++) {
		int temp = WM[i][h] + WM(h+1,j);
		if (temp < WMijp)
			WMijp = temp;
	}
	WMij = Vij + auPen(rnai, rnaj) + b;
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(RNA[i + 1], rnaj) + b + c;
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(rnai, RNA[j - 1]) + b + c;
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1] + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1], RNA[j - 1]) + b + 2* c ;

	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));
	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;
	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];
	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];
	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WMij = MIN(WMijp, WMij);
	WM[i][j] = WMij;
	WM(i,j) = WMij;
	/* WM finaliza*/
	return;
}
/* Función para calcular WM, V y WM en un punto (i,j). El calculo de V requiere VBI y WM, además el calculo de WM requiere V */
void calcVMVWM(int i, int j) {
	int a = multConst[0] /*penalidad para los multiloops*/, 
	b = multConst[2]/*penalidad por rama para los multiloops*/, 
	c = multConst[1]/* Penalidad por base singular para los multiloops*/, a1, h, es;
	int aupen;
	int WMijp, WMidjd, WMidj, WMijd, WMij;
	int VMij, VMijd, VMidj, VMidjd, A_temp;
	int WMip1hm1 /* Valor de WM con i+1 y j-1 */,
	WMip2hm1 /* Valor de WM con i+2 y h-1*/,
	WMhjm1/* Valor de WM en h con j-1*/, WMhjm2/* Valor de WM en h con j-2*/,
	WMhp1j /*Valor de WM en h+1 con j*/;
	int rnai, rnaj;
	int tmp1, tmp2;
	rnai = RNA[i];
	rnaj = RNA[j];
	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;
	/* VM y WM inician */
	aupen = auPen(rnai, rnaj); /* penalidad para el par de bases (i,j) que sea AU o no GC */
	VMij = INFINITY_;
	VMijd = INFINITY_;
	VMidj = INFINITY_;
	VMidjd = INFINITY_;
	WMijp = WM[i][i + 4] + WM(i+5,j);
	a1 = WM[i][i + 5] + WM(i+6,j);
	if (a1 <= WMijp)
		WMijp = a1;
	/* Aquí se hace el calculo de VM y WM en (i,j) concurrentemente 
	. Este for calcula los valores de VMij, VMidj, VMidjd usando WMijp.
	El valor de WMijp es necesitado para calcular el valor de WM en (i,j). De cualquier forma, debiera ser notado
	que el valor final de WM[i][j] requiere de V en (i,j) y es calculado al final 
	*/
	/* Hay 4 posibilidades para el multiloop cerrando el par de bases relacionado con la inclusion de energias colgantes */
	/* 1)Incluyendo la energia colgante de la base i+1 y de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMidjd
	   2)Incluyendo la energia colgante de la base i+1 y no la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMidj
	   3)No Incluyendo la energia colgante de la base i+1 y si la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMijd	 
	   4)No Incluyendo la energia colgante de la base i+1 y tampoco la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMij */
	for (h = i + 6; h <= j - 5; h++) {
		a1 = WM[i][h];
		WMip1hm1 = WM[i + 1][h - 1];
		WMip2hm1 = WM[i + 2][h - 1];
#if 0
		WMhjm1 = WM[h][j-1];
		WMhjm2 = WM[h][j-2];
		WMhp1j = WM[h+1][j];
#else
		WMhjm1 = WM(h,j-1); /* El preprocesador convertirá esto en WM[j-1][h]. De acuerdo al algoritmo
		  la expresion debiera ser WM[h][j-1], en ese caso como el valor de h cambia, esto permitirá acceder
		  a elementos de una columna de la matriz WM. Para mejorar tiempo de ejecución la matriz WM se hace
		  simétrica, demanera que WM[h][j-1] = WM[j-1][h] y entonces tendra accesos en una fila  */ 
		WMhjm2 = WM(h,j-2); 
		WMhp1j = WM(h+1,j); 
#endif
		/* WM inicia */
		a1 += WMhp1j;
		if (a1 <= WMijp)
			WMijp = a1;
		/* WM finaliza */
		/* Calculo de las 4 opciones para VM*/
		A_temp = WMip1hm1 + WMhjm1;
		if ((A_temp <= VMij))
			VMij = A_temp;
		A_temp = WMip2hm1 + WMhjm1;
		if (A_temp <= VMidj && constraints[i + 1] <= 0)
			VMidj = A_temp;
		A_temp = WMip1hm1 + WMhjm2;
		if (A_temp <= VMijd && constraints[j - 1] <= 0)
			VMijd = A_temp;
		A_temp = WMip2hm1 + WMhjm2;
		if (A_temp <= VMidjd && constraints[i + 1] <= 0 && constraints[j - 1] <= 0)
			VMidjd = A_temp;
	}

#if 0
	VMidj += dangle[rnai][rnaj][RNA[i+1]][0];
	VMijd += dangle[rnai][rnaj][RNA[j-1]][1];
	VMidjd += dangle[rnai][rnaj][RNA[i+1]][0] + dangle[rnai][rnaj][RNA[j-1]][1];
#else
	tmp1 = dangle[rnai][rnaj][RNA[i + 1]][0]; /* energia colgante del par de bases (i,j) con una base singular i+1 en el final 5' */ 
	tmp2 = dangle[rnai][rnaj][RNA[j - 1]][1]; /* energia colgante del par de bases (i,j) con una base singular i+1 en el final 5' */
	VMidj += (tmp1 + c);
	VMidjd += (tmp1 + c);
	VMijd += (tmp2 + c);
	VMidjd += (tmp2 + c);
#endif
	VMij = MIN(MIN(VMij, VMidj), MIN(VMijd, VMidjd));
	VMij = VMij + b + a + aupen;
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		VMij = INFINITY_;
	VM[i][j] = VMij;
	/* VM finaliza */
	/* V inicia */
	es = eS(i, j);
	if (es == 0) {
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];
	int Vij;
	Vij = MIN(MIN(eH(i, j), es), MIN(VBI[i][j], VMij));
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 	&& constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		Vij = INFINITY_;
	V[indx[i] + j] = Vij;
	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		//Se revisa si i+1, j-1 estan pareadas
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			//si no lo están se revisa por i-1,j+1
			//los pares de bases isolados se ven mas adelante 
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			if (Vijl > INFINITY_ - SMALLINFTY_)
				//par de bases isoladas encontrada, se setea la energía a infinito
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		printf("ERROR: En calcVBIVMVWM: i: %5d  j: %5d\n",i,j);
#endif
	/* V finaliza */
	/* WM inicia */
	//Se trabaja sobre los WMs
	WMij = Vij + auPen(rnai, rnaj) + b;
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(
				RNA[i + 1], rnaj) + b + c;
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(
				rnai, RNA[j - 1]) + b + c;
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1] + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1], RNA[j - 1]) + b + 2* c ;

	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));
	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;
	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];
	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];
	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WM[i][j] = MIN(WMijp, WMij);
	WM(i,j) = WM[i][j]; 
	/* WM finaliza */
	return;
}

/* Calculo de VBI, VM V y WM. El calculo de V en (i,j) requiere VBI y VM. Además, el calculo del valor final
   de WM(i,j) requiere V en (i,j)  */
void calcVBIVMVWM(int i, int j) {
	int ip, jp, temp, VBIij, thres1;
	int
	a = multConst[0] /* a es una penalidad para multiloops  */,
	b = multConst[2]/* b es una penalidad para las ramas de multiloops, una por rama  */,
	c = multConst[1] /* c es una penalidad para nucleotidos monocaterianos en los multiloops  */,
	a1, h, es;
	int aupen;
	int WMijp, WMidjd, WMidj, WMijd, WMij;
	int VMij, VMijd, VMidj, VMidjd, A_temp;
	int WMip1hm1 /* Valor de WM en i+1 y h-1 */,
	WMip2hm1 /* Valor de WM en i+2 y j-1*/,
	WMhjm1 /* Valor de WM en h y j-1*/,
	WMhjm2 /* Valor de WM en h y j-2*/, WMhp1j /*Valor de WM en h+1 y j*/;
	int rnai, rnaj;
	int tmp1, tmp2;
	rnai = RNA[i];
	rnaj = RNA[j];
	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;
	/* VBI inicia */
	VBIij = INFINITY_;
	int ifinal, jfinal;
	/* ip=i+1, jp=j-1 cierran un stack, por lo que el limite de jp se define como j-1 en el siguiente loop  */
	ip = i + 1;
	thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); /* El tamaño minimo de un loop de horquilla es 3. Por ello, jp inicia de ip+4  */
	for (jp = thres1; jp <= j - 2; jp++) {
		if (chPair(RNA[ip], RNA[jp])) {
			if (checkSS(i, ip) || checkSS(jp, j))
				continue;
			temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; /* Energia de un loop interno cerrado por (i,j)  
			 y (ip,jp) + la energia de la estructura cerrada por (ip,jp)  */ 
			if (VBIij > temp) {
				VBIij = temp;
				ifinal = ip;
				jfinal = jp;
			}
		}
	}
	for (ip = i + 2; ip <= i + MAXLOOP + 1; ip++) {
		thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4);
		for (jp = thres1; jp <= j - 1; jp++) {
			if (chPair(RNA[ip], RNA[jp])) {
				if (checkSS(i, ip) || checkSS(jp, j))
					continue;
				temp = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > temp) {
					VBIij = temp;
					ifinal = ip;
					jfinal = jp;
				}
			}
		}
	}
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		VBIij = INFINITY_;
	VBI[i][j] = VBIij;
	/* VBI finaliza */
	/* VM and WM inician */
	aupen = auPen(rnai, rnaj);
	/* Hay 4 posibilidades para el multiloop cerrando el par de bases relacionado con la inclusion de energias colgantes */
	/* 1)Incluyendo la energia colgante de la base i+1 y de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMidjd
	   2)Incluyendo la energia colgante de la base i+1 y no la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMidj
	   3)No Incluyendo la energia colgante de la base i+1 y si la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMijd	 
	   4)No Incluyendo la energia colgante de la base i+1 y tampoco la de la base j-1 con el par de bases (i,j) cerrando el multiloop - WMij */

	VMij = INFINITY_;
	VMijd = INFINITY_;
	VMidj = INFINITY_;
	VMidjd = INFINITY_;
	/* Calculos combinados de WM(i,j) y VM(i,j)  */
	WMijp = WM[i][i + 4] + WM(i+5,j);
	a1 = WM[i][i + 5] + WM(i+6,j);
	if (a1 <= WMijp)
		WMijp = a1;
	for (h = i + 6; h <= j - 5; h++) {
		a1 = WM[i][h];
		WMip1hm1 = WM[i + 1][h - 1];
		WMip2hm1 = WM[i + 2][h - 1];
#if 0
		WMhjm1 = WM[h][j-1];
		WMhjm2 = WM[h][j-2];
		WMhp1j = WM[h+1][j];
#else
		WMhjm1 = WM(h,j-1);
		WMhjm2 = WM(h,j-2);
		WMhp1j = WM(h+1,j);
#endif
		/* WM inicia */
		a1 += WMhp1j;
		if (a1 <= WMijp)
			WMijp = a1;
		A_temp = WMip1hm1 + WMhjm1;
		if ((A_temp <= VMij))
			VMij = A_temp;
		A_temp = WMip2hm1 + WMhjm1;
		if (A_temp <= VMidj && constraints[i + 1] <= 0)
			VMidj = A_temp;
		A_temp = WMip1hm1 + WMhjm2;
		if (A_temp <= VMijd && constraints[j - 1] <= 0)
			VMijd = A_temp;
		A_temp = WMip2hm1 + WMhjm2;
		if (A_temp <= VMidjd && constraints[i + 1] <= 0 && constraints[j - 1] <= 0)
			VMidjd = A_temp;
	}

#if 0
	VMidj += dangle[rnai][rnaj][RNA[i+1]][0];
	VMijd += dangle[rnai][rnaj][RNA[j-1]][1];
	VMidjd += dangle[rnai][rnaj][RNA[i+1]][0] + dangle[rnai][rnaj][RNA[j-1]][1];
#else
	tmp1 = dangle[rnai][rnaj][RNA[i + 1]][0];
	tmp2 = dangle[rnai][rnaj][RNA[j - 1]][1];

	VMidj += tmp1;
	VMidjd += tmp1;
	VMijd += tmp2;
	VMidjd += tmp2;
#endif
	VMij = MIN(MIN(VMij, VMidj), MIN(VMijd, VMidjd));
	VMij = VMij + b + a;
	VMij += aupen;
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		VMij = INFINITY_;
	VM[i][j] = VMij;
	es = eS(i, j); /* Energía de un stack cerrado con (i,j) y (i+1, j-1)*/
	if (es == 0) { 
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];
	int Vij;
	Vij = MIN(MIN(eH(i, j), es), MIN(VBI[i][j], VMij));
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0 && constraints[j] != i) || constraints[i] == -1 || constraints[j] == -1)
		Vij = INFINITY_;
	V[indx[i] + j] = Vij;
	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			if (Vijl > INFINITY_ - SMALLINFTY_)
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		printf("ERROR: En calcVBIVMVWM: i: %5d  j: %5d\n",i,j);
#endif
	/* V finaliza */
	/* WM inicia */
	/* No hay base colgante en ningun lado */
	WMij = V[indx[i] + j] + aupen + b;
	/* Base colgante i en el final 3'del par de bases (i+1,j) */
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(RNA[i + 1], rnaj) + b + c;
	/* Base colgante j en el final 5' del par de bases (i,j-1)  */
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(rnai, RNA[j - 1]) + b + c;
	/* Base colgante ien el final 3' y base j en el final 5' del par de bases (i+1, j-1)  */  
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1] + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1], RNA[j - 1]) + b + 2* c ;
	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));
	/* El termino WM[i+1][j] se encarga cuando la base i no se parea o no esta cumpliendo el rol de calcular
	   la energía colgante y se debe agregar la penalidad c para que la base i se mantenga singular */
	/* El termino WM[i][j-1] se encarga cuando la base j no se parea o no esta cumpliendo el rol de calcular
	   la energía colgante y se debe agregar la penalidad c para que la base j se mantenga singular */
	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;
	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];
	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];
	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WM[i][j] = MIN(WMijp, WMij);
	WM(i,j) = WM[i][j]; 
	/* WM finaliza */
	return;
}
///
/* Función para calcular el valor de W[j] */
void calcW(int j) {
	int i;
	int Wj, Widjd /* Base colgante en amos lados */,
		Wijd/* Base colgante en el j-esimo lado.*/,
		Widj/* Base colgante en el i-esimo lado */,
		Wij/* Sin bases colgantes en ningun lado */, Wim1 /* Valor de W en (i-1). Se setea en cero si es positivo*/;
	int rnai, rnaj;
	int must_branch = 0, besti = 0;
	Wj = INFINITY_;
	rnaj = RNA[j];
	for (i = 1; i < j - 3; i++) {
		Wij = Widjd = Wijd = Widj = INFINITY_;
# if 0
		Wim1=W[i-1];
#endif
#if 1
		Wim1 = MIN(0, W[i - 1]);/* Si W[i-1]>=0, esto significa que hay una rama contenida en secuencia
		que va de 1 hasta i-1. En cualquier otro caso W[i-1] será infinito. Aquí Wim1 está definido de esta
		manera, para hacer que la enrgía de una secuencia sin plegar sea infinito*/ 
		
#endif
		//Wim1 = W[i - 1];
		rnai = RNA[i];
		/* Se calcula la energía sin bases colgantes  */
		Wij = V[indx[i] + j] + auPen(rnai, rnaj) + Wim1;
		/* Se cuelga en ambos lados del par de bases (i+1,j-1). Se agrega la energía correspondiente  */
		if (constraints[i] <= 0 && constraints[j] <= 0)
			Widjd = V[indx[i + 1] + j - 1] + auPen(RNA[i + 1], RNA[j - 1]) + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j- 1]][RNA[i + 1]][rnaj][0] + Wim1;
		/* Base singular j colgando en el final 5' del par de bases (i,j-1) */
		if (constraints[j] <= 0)
			Wijd = V[indx[i] + j - 1] + auPen(rnai, RNA[j - 1]) + dangle[RNA[j- 1]][rnai][rnaj][0] + Wim1;
		/* Base singular i colgando en el final 3' del par de bases (i+1,j) */
		if (constraints[i] <= 0)
			Widj = V[indx[i + 1] + j] + auPen(RNA[i + 1], rnaj) + dangle[rnaj][RNA[i + 1]][rnai][1] + Wim1;
		int tmpWj = Wj;
		Wj = MIN(MIN(MIN(Wij, Widjd), MIN(Wijd, Widj)), Wj); /* Se toma el mínimo */
		if (tmpWj != Wj) {
			must_branch = 0;
			besti = i;
		}
		if (Wj < INFINITY_) {
			if (Wj == Wij) {
				if (constraints[i] == j) {
					must_branch = 1;
				}
			} else if (Wj == Widjd) {
				if (constraints[i + 1] == j - 1) {
					must_branch = 1;
				}
			} else if (Wj == Wijd) {
				if (constraints[i] == j - 1) {
					must_branch = 1;
				}
			} else {
				if (constraints[i + 1] == j) {
					must_branch = 1;
				}
			}
		}
	}
	/* Si la j-esima base no esta contribuyendo en el calculo de energía de W[j]*/
	if (!must_branch) {
		if (Wj > W[j - 1])
			Wj = W[j - 1];
	}
	W[j] = Wj;
	return;
}

/* Se calcula el valor de energía del loop interno con (i,j) como par de bases cerrando y (ip,jp) como par de
bases encerrado*/
int eL(int i, int j, int ip, int jp) {
	int energy;
	int size1, size2, size;
	int loginc; 
	int lopsided; /* Define la asimetría de un loop interior */
	energy = INFINITY_;
	loginc = 0;
	size1 = ip - i - 1;
	size2 = j - jp - 1;
	size = size1 + size2;
	if (size1 == 0 || size2 == 0) {
		if (size > 30) {
			/*  No depende de i,j,ip o jp */
			loginc = (int) floor(prelog * log((double) size / 30.0));
			energy = bulge[30] + eparam[2] + loginc + auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]);
		} else if (size <= 30 && size != 1) {
			/*  No depende de i,j,ip o jp */
			energy = bulge[size] + eparam[2];
			energy += auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]);
		} else if (size == 1) {
			energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[ip], RNA[jp])] + bulge[size] + eparam[2]; /* mans */
		}
	} else {
		/* Loop interno */
		lopsided = abs(size1 - size2);
		if (size > 30) {
			loginc = (int) floor(prelog * log((double) size / 30.0));
			if (!((size1 == 1 || size2 == 1) && gail)) { /* Loop interno normal con tamaño >30*/
				energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j- 1])] + tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp+ 1], RNA[ip - 1])] + inter[30] + loginc + eparam[3]+ MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1,size2))]));
			} else { 
				energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)] + inter[30] + loginc + eparam[3] + MIN(maxpen, (lopsided	* poppen[MIN(2, MIN(size1, size2))]));
			}
		}
		/* si el tamaño no es mayor a 30 se generan muchos casos */
		else if (size1 == 2 && size2 == 2) {
			/* loop interno de 2x2 */
			energy = iloop22[RNA[i]][RNA[ip]][RNA[j]][RNA[jp]][RNA[i + 1]][RNA[i+ 2]][RNA[j - 1]][RNA[j - 2]];
		} else if (size1 == 1 && size2 == 2) {
			energy = iloop21[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]][RNA[j - 2]][RNA[ip]][RNA[jp]];
		} else if (size1 == 2 && size2 == 1) {
			/* Loop interno de 1x2 */
			energy = iloop21[RNA[jp]][RNA[ip]][RNA[j - 1]][RNA[i + 2]][RNA[i+1]][RNA[j]][RNA[i]];
		} else if (size == 2) {
			/* Loop interno de 1x1*/
			energy = iloop11[RNA[i]][RNA[i + 1]][RNA[ip]][RNA[j]][RNA[j - 1]][RNA[jp]];
		} else if ((size1 == 1 || size2 == 1) && gail) { 
			energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)] + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)] + inter[size] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
		} else { 
			energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[size] + loginc + eparam[3] + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1, size2))]));
		}
	}
	return energy;
}

/* Funcion usada para calcular la energia de un loop de horquilla entre i & j*/
int eH(int i, int j) {
	/*  Loop de horquilla para todas las bases entre i y j*/
	/*  size es el tamaño del loop, energy es el resultado, loginc es para la extrapolacion de loops mayores a 30*/
	int size;
	int loginc;
	int energy = INFINITY_;
	int key, index, count, tlink, kmult;
	size = j - i - 1; //size es el numero de bases en el loop cuando el par que cierra es excluido 
	//se verifica so la region monocateriana esta permitida con las constantes dadas
	if (checkSS(i, j) || (constraints[i] > 0 && constraints[i] != j)
			|| (constraints[j] > 0 && constraints[j] != i))
		return energy;
	/* Se ve en hairpin el valor correcto de los 30 presentes*/
	if (size > 30) {
		loginc = (int) ((prelog) * log(((double) size) / 30.0));
		energy = hairpin[30] + loginc + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4]; /* Penalidad de tamaño + energia de stacking incompatible*/
	}
	else if (size <= 30 && size > 4) {
		energy = hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4]; /* Penalidad de tamaño + energia de stacking incompatible*/
	}
	else if (size == 4) {
		/*  tetraloop */
		key = 0;
		tlink = 0;
		for (index = 0; index < 6; ++index) {
			switch (RNA[i + index]) {
			case BASE_A:
				kmult = 1;
				break;
			case BASE_C:
				kmult = 2;
				break;
			case BASE_G:
				kmult = 3;
				break;
			case BASE_U:
				kmult = 4;
				break;
			default:
				kmult = 0;
				fprintf(stderr, "ERROR: in tetraloop calculation\n");
			}
			key += kmult * (int) pow(10.0, 5 - index);
		}
		/* Si la secuencia está en el tetraloop se usa este valor*/
		for (count = 1; count < numoftloops && tlink == 0; ++count) {
			if (key == tloop[count][0]) {
				tlink = tloop[count][1];
			}
		}
		energy = tlink + hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4];
	}
	else if (size == 3) {
		energy = hairpin[size];
		energy += auPen(RNA[i], RNA[j]);
	}
	else if (size < 3 && size != 0) {
		/*  Sin incopatibilidad terminal */
		energy = hairpin[size] + eparam[4];
		if ((RNA[i] == BASE_A && RNA[j] == BASE_U) || (RNA[i] == BASE_U
				&& RNA[j] == BASE_A)) {
			energy += 6; 
			/* Loops de horquilla de tamaño 3 no son permitidos*/
		}
	} else if (size == 0)
		return INFINITY_;
	/*  i-2 = i-1 = i = G, and j = U; i < j */
	if (i > 2) {
		if (RNA[i - 2] == BASE_G && RNA[i - 1] == BASE_G && RNA[i] == BASE_G
				&& RNA[j] == BASE_U) {
			energy += gubonus;
		}
	}
	tlink = 1;
	for (index = 1; (index <= size) && (tlink == 1); ++index) {
		if (RNA[i + index] != BASE_C)
			tlink = 0;
	}
	if (tlink == 1) {
		if (size == 3) {
			energy += c3;
		} else {
			energy += cint + size * cslope;
		}
	}
	return energy;
}

/* Función usada para calcular la energía de pares de stack (i,j) & (i+1,j-1)*/
int eS(int i, int j) {
	int energy;
	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i))
		return INFINITY_;
	energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])] + eparam[1];
	return energy;
}
