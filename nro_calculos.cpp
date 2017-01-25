/*
	Mini-Algoritmo que permite saber a cuantas secuencias se les va a calcular el ZScore en el programa principal (main.cpp)
 */

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

using namespace std;

int main(){
	int secuencia=0,ventana=0,corredor=0,contador=0;
	cout << "Ingrese el largo de la secuencia: ";
	cin >> secuencia;
	cout << "Ingrese el largo de la ventana: ";
	cin >> ventana;
	while (corredor+ventana<=secuencia)
	{
		contador++;
		corredor = corredor + ventana/2;
	}
	cout << "La cantidad de segmentos a analizar es: "<<contador<<endl;
	return 0;	
}
