#include "klu_internal.h"
#include <string.h>
#define PRECISION 20

void dumpKPerm(int* Q, int* P, int n, int counter)
{
    int i;
    char strQ[32];
    char strQi[32];
    char strP[32];
    char strPi[32];
    char counterstring[32];
    sprintf(counterstring, "%d", counter);

    strcpy(strQ, "KLU_Q");
    strcpy(strP, "KLU_P");
    strcat(strQ, counterstring);
    strcat(strP, counterstring);
    strcat(strQ, ".txt");
    strcat(strP, ".txt");

    /* dump permutations into file */
    FILE* fQ = fopen(strQ, "w");
    FILE* fP = fopen(strP, "w");

    /* dump Q */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fQ, "%d, ", Q[i]);
    }
    fprintf(fQ, "%d\n", Q[n-1]);

    /* dump P */
    for (i = 0 ; i < n-1; i++)
    {
        fprintf(fP, "%d, ", P[i]);
    }
    fprintf(fP, "%d\n", P[n-1]);

    fclose(fQ);
    fclose(fP);
}

void printKMTX(FILE* f, double* values, int* indices, int* pointers, int n, int nz)
{
    int i,j;
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%d %d %d\n", n, n, nz);
    for(i = 0; i < n ; i++)
    {
        for(j = pointers[i] ; j < pointers[i]+1 ; j++)
        {
            fprintf(f, "%d %d %.*e\n", indices[j]+1, i+1, PRECISION, values[j]);
        }
    }
}

void dumpKLU(double *Lx, int *Li, int *Lp, 
            double *Ux, int *Ui, int *Up, 
            double *Fx, int *Fi, int *Fp, 
            int lnz, int unz, int n, 
            int nzoff, int counter)
{
    char strL[32];
    char strU[32];
    char strF[32];
    char counterstring[32];
    sprintf(counterstring, "%d", counter);
    strcpy(strL, "KLU_L");
    strcpy(strU, "KLU_U");
    strcpy(strF, "KLU_F");
    strcat(strL, counterstring);
    strcat(strU, counterstring);
    strcat(strF, counterstring);
    strcat(strL, ".mtx");
    strcat(strU, ".mtx");
    strcat(strF, ".mtx");

    FILE *l, *u, *g;
    g = fopen(strF, "w");
    l = fopen(strL, "w");
    u = fopen(strU, "w");

    printKMTX(g, Fx, Fi, Fp, n, nzoff);
    printKMTX(l, Lx, Li, Lp, n, lnz);
    printKMTX(u, Ux, Ui, Up, n, unz);

    fclose(g);
    fclose(l);
    fclose(u);
}

void dumpKPath(int* path, int* bpath, int n, int nb, int counter)
{
        int i;
        char counterstring[32];
        sprintf(counterstring, "%d", counter);
        char strbpath[32];
        char strpath[32];
        strcpy(strbpath, "KLU_bpath");
        strcpy(strpath, "KLU_path");
        strcat(strbpath, counterstring);
        strcat(strpath, counterstring);
        strcat(strbpath, ".txt");
        strcat(strpath, ".txt");
        FILE* fbpath = fopen(strbpath, "w");
        FILE* fpath = fopen(strpath, "w");

        /* dump path into file */
        for (i = 0; i < nb - 1; i++)
        {
            fprintf(fbpath, "%d, ", bpath[i]);
        }
        fprintf(fbpath, "%d\n", bpath[nb-1]);

        for (i = 0; i < n - 1; i++)
        {
            fprintf(fpath, "%d, ", path[i]);
        }
        fprintf(fpath, "%d\n", path[n-1]);

        fclose(fpath);
        fclose(fbpath);
}

void dumpKAll(double *Lx, 
            int *Li, 
            int *Lp,
            double *Ux, 
            int *Ui, 
            int *Up, 
            double *Fx, 
            int *Fi, 
            int *Fp, 
            int *P,
            int *Q,
            int *path,
            int *bpath,
            int lnz,
            int unz,
            int n,
            int nzoff,
            int nb
        )
{
    static int counter = 0;
    dumpKPerm(Q, P, n, counter);
    dumpKLU(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, lnz, unz, n, nzoff, counter);
    dumpKPath(path, bpath, n, nb, counter);
    counter++;
}