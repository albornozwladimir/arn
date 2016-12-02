
#ifndef _MAIN_C_H
#define _MAIN_C_H

#define DYNALLOC

#ifdef DYNALLOC
extern int LENGTH;
extern unsigned char *RNA1; /* [LENGTH] */
extern unsigned char *RNA; /* [LENGTH] */
extern int *structure; /* [LENGTH] */
extern int *V; /* [(LENGTH-1)*(LENGTH)/2 + 1] */
extern int *W; /* [LENGTH] */
extern int **VBI; /* [LENGTH][LENGTH] */
extern int **VM; /* [LENGTH][LENGTH] */
extern int **WM; /* [LENGTH][LENGTH] */
extern int *indx; /* [LENGTH] */
extern int *constraints;
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

#endif
