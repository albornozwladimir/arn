#ifndef _CONSTANTS_H
#define _CONSTANTS_H

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

#endif
