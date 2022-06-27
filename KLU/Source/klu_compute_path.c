/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/*
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"
#include <string.h>

void dumpPerm(int* Q, int* P, int n, int counter)
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
    // strcat(strQi, counterstring);
    // strcat(strPi, counterstring);
    // strcpy(strQi, "KLU_Qi");
    // strcpy(strPi, "KLU_Pi");
    // strcat(strQi, ".txt");
    // strcat(strPi, ".txt");

    /* dump permutations into file */
    FILE* fQ = fopen(strQ, "w");
    FILE* fP = fopen(strP, "w");
    // FILE* fQi = fopen(strQi, "w");
    // FILE* fPi = fopen(strPi, "w");

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
    // fclose(fQi);
    // fclose(fPi);
}

void dumpLU(double *Lx, int *Li, int *Lp, double *Ux, int *Ui, int *Up, double *Fx, int *Fi, int *Fp, int lnz, int unz,
            int n, int nzoff, int counter)
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
    strcat(strL, ".csc");
    strcat(strU, ".csc");
    strcat(strF, ".csc");

    FILE *l, *u, *g;
    g = fopen(strF, "w");
    l = fopen(strL, "w");
    u = fopen(strU, "w");
    int i;

    /* Print off-diagonal blocks in csr format */
    for (i = 0; i < nzoff - 1; i++)
    {
        fprintf(g, "%.17lf, ", Fx[i]);
    }
    fprintf(g, "%lf\n", Fx[nzoff-1]);
    for (i = 0; i < nzoff-1; i++)
    {
        fprintf(g, "%d, ", Fi[i]);
    }
    fprintf(g, "%d\n", Fi[nzoff-1]);
    for (i = 0; i < n ; i++)
    {
        fprintf(g, "%d, ", Fp[i]);
    }
    fprintf(g, "%d\n", Fp[n]);

    /* Print L-matrix in csr format */
    for (i = 0; i < lnz-1; i++)
    {
        fprintf(l, "%lf, ", Lx[i]);
    }
    fprintf(l, "%.17lf\n", Lx[lnz-1]);
    for (i = 0; i < lnz-1; i++)
    {
        fprintf(l, "%d, ", Li[i]);
    }
    fprintf(l, "%d\n", Li[lnz-1]);
    for (i = 0; i < n; i++)
    {
        fprintf(l, "%d, ", Lp[i]);
    }
    fprintf(l, "%d\n", Lp[n]);

    /* Print U-matrix in csr format */
    for (i = 0; i < unz-1; i++)
    {
        fprintf(u, "%.17lf, ", Ux[i]);
    }
    fprintf(u, "%lf\n", Ux[unz-1]);
    for (i = 0; i < unz-1; i++)
    {
        fprintf(u, "%d, ", Ui[i]);
    }
    fprintf(u, "%d\n", Ui[unz-1]);
    for (i = 0; i < n; i++)
    {
        fprintf(u, "%d, ", Up[i]);
    }
    fprintf(u, "%d\n", Up[n]);

    fclose(l);
    fclose(u);
    fclose(g);
}

void dumpPath(int* path, int* bpath, int n, int nb, int counter)
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

void dumpAll(double *Lx, 
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
    dumpPerm(Q, P, n, counter);
    dumpLU(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, lnz, unz, n, nzoff, counter);
    dumpPath(path, bpath, n, nb, counter);
    counter++;
}

/*
 * Main function. Expects a factorized matrix and "changeVector", which contains columns of A that change
 * e.g. changeVector = {3, 5}, if columns 3 and 5 contain varying entries. Needs to be permuted, since
 * L*U = P*A*Q
 */
int KLU_compute_path(KLU_symbolic *Symbolic, KLU_numeric *Numeric, KLU_common *Common, Int *changeVector, Int changeLen)
{
    /* Declarations */
    /* LU data */
    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int nb = Symbolic->nblocks;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double *Rs;
    // int* Pinv = Numeric->Pinv;

    int RET;

    /* indices and temporary variables */
    int i, k, j, ent;
    int pivot;
    int col, nextcol;
    int u_closest, l_closest;

    // blocks
    Int k2, k1, nk;
    // TODO: save sizeof(...) statically and not call for each alloc
    Int *Qi = calloc(n, sizeof(Int));
    Int *changeVector_permuted = calloc(changeLen, sizeof(Int));

    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->bpath)
    {
        KLU_free(Numeric->bpath, Numeric->nblocks, sizeof(int), Common);
    }

    Numeric->path = KLU_malloc(n, sizeof(int), Common);
    Numeric->bpath = KLU_malloc(nb, sizeof(int), Common);

    for (i = 0; i < n; i++)
    {
        Numeric->path[i] = 0;
    }
    for (i = 0; i < nb; i++)
    {
        Numeric->bpath[i] = 0;
    }

    // TODO: solve smarter, no more klu_extracts.
    Lp = calloc(n + 1, sizeof(int));
    Up = calloc(n + 1, sizeof(int));
    Fp = calloc(n + 1, sizeof(int));
    Lx = calloc(lnz, sizeof(double));
    Ux = calloc(unz, sizeof(double));
    Fx = calloc(nzoff, sizeof(double));
    Li = calloc(lnz, sizeof(int));
    Ui = calloc(unz, sizeof(int));
    Fi = calloc(nzoff, sizeof(int));
    P = calloc(n, sizeof(int));
    Q = calloc(n, sizeof(int));
    Rs = calloc(n, sizeof(double));
    R = calloc(nb + 1, sizeof(int));

    // first, get LU decomposition
    // sloppy implementation, as there might be a smarter way to do this
    RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);
    //saveLU(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, lnz, unz, n, nzoff);
    if (RET != (TRUE))
    {
        return (FALSE);
    }

    // second, invert permutation vector
    // Q gives "oldcol", we need "newcol"
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    // third, apply permutation on changeVector
    for (i = 0; i < changeLen; i++)
    {
        changeVector_permuted[i] = Qi[changeVector[i]];
    }

    for(i=0 ; i < changeLen ; i++)
    {
        printf("%d, ", changeVector_permuted[i]);
    }

    // fourth, sort permuted vector
    // in "full partial refactorisation", only find minimum value
    // sloppy selectionsort implementation for first design
    // for (i = 0 ; i < changeLen ; i++)
    // {
    //     k = i;
    //     pivot = changeVector_permuted [i];
    //     for (j = i + 1 ; j < changeLen ; j++)
    //     {
    //         if (changeVector_permuted [j] < pivot)
    //         {
    //             pivot = changeVector_permuted [j];
    //             k = j;
    //         }
    //     }
    //     if (k != i)
    //     {
    //         // found smaller
    //         // switch positions
    //         int tmp = changeVector_permuted [i];
    //         changeVector_permuted [i] = pivot;
    //         changeVector_permuted [k] = tmp;
    //     }
    // }

    // step three and four can / should be done externally

    if (Common->btf == FALSE)
    {
        Numeric->bpath[0] = 1;
        /* no blocks */
        for (i = 0; i < changeLen; i++)
        {
            pivot = changeVector_permuted[i];
            if (Numeric->path[pivot] == 1)
            {
                continue;
            }
            Numeric->path[pivot] = 1;
            while (pivot < n)
            {

                u_closest = n + 1;
                l_closest = n + 1;

                // find closest off-diagonal entry in L
                col = Lp[pivot];
                nextcol = Lp[pivot + 1];

                if (nextcol - col == 1)
                {
                    // only one entry in column => diagonal entry
                    l_closest = n + 1;
                }
                else
                {
                    // indices are not sorted!!!! TODO. Use klu_sort maybe?

                    // find closest off-diagonal entry in L, "look down" in pivot-th column
                    for (ent = col + 1; ent < nextcol; ent++)
                    {
                        if (l_closest > Li[ent])
                        {
                            l_closest = Li[ent];
                        }
                    }
                }
                if (l_closest - pivot == 1)
                {
                    // l_closest = pivot + 1, there can't be a closer off-diagonal row entry. save computation time.
                    goto minimum;
                }

                // find closest off-diagonal entry in U
                /* u saved in csc:
                   Ux = [x, x, x, x, x] entries saved column-wise
                   Ui = [0, 2, 1, 4, 3] row position of those entries in columns
                   Up = [0, 1, 2, 3, 4, 5] column pointers
                   to find closest off-diagonal of row k:
                   find all Ui[j] == pivot
                   find Up[k] <= j < Up[k+1], k > pivot
                   ALTERNATIVELY
                   Transpose U
                   do the same as for L
                 */
                for (j = 0; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        // check if Ui[j] is right to pivot
                        for (k = 0; k < n + 1; k++)
                        {
                            if (Up[k] <= j && j < Up[k + 1])
                            {
                                if (k > pivot)
                                {
                                    // found closest off-diagonal
                                    u_closest = k;
                                    // don't hate me
                                    goto minimum;
                                    // alternatively: j=unz+1
                                }
                                break;
                            }
                        }
                    }
                }
            minimum:
                pivot = MIN(l_closest, u_closest);
                // check if pivot is either already in path or n+1 (no more off-diag values)
                if (Numeric->path[pivot] == 1 || pivot == n + 1)
                {
                    break;
                }
                Numeric->path[pivot] = 1;
            }
        }
    }
    else
    {
        // fifth, compute factorization path
        for (i = 0; i < changeLen; i++)
        {
            // get next changing column
            pivot = changeVector_permuted[i];

            // check if it was already computed
            if (Numeric->path[pivot] == 1)
            {
                // already computed pivot
                // do nothing, go to next pivot
                continue;
            }

            // set first value of singleton path
            Numeric->path[pivot] = 1;

            // find block of pivot
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */
                    Numeric->bpath[k] = 1;
                    break;
                }
            }

            if (nk == 1)
            {
                // 1x1-block, its pivot already in path
                continue;
            }

            // propagate until end
            // in blocks, pivot < n_block[k]
            while (pivot < k2)
            {
                u_closest = k2;
                l_closest = k2;

                // find closest off-diagonal entry in L
                col = Lp[pivot];
                nextcol = Lp[pivot + 1];

                if (nextcol - col == 1)
                {
                    // only one entry in column => diagonal entry
                    l_closest = k2;
                }
                else
                {
                    // indices are not sorted!!!! TODO. Use klu_sort maybe?

                    // find closest off-diagonal entry in L, "look down" in pivot-th column
                    for (ent = col + 1; ent < nextcol; ent++)
                    {
                        if (l_closest > Li[ent])
                        {
                            l_closest = Li[ent];
                        }
                    }
                }
                if (l_closest - pivot == 1)
                {
                    // l_closest = pivot + 1, there can't be a closer off-diagonal row entry. save computation time.
                    goto minimum_btf;
                }

                // find closest off-diagonal entry in U
                /* u saved in csc:
                   Ux = [x, x, x, x, x] entries saved column-wise
                   Ui = [0, 2, 1, 4, 3] row position of those entries in columns
                   Up = [0, 1, 2, 3, 4, 5] column pointers
                   to find closest off-diagonal of row k:
                   find all Ui[j] == pivot
                   find Up[k] <= j < Up[k+1], k > pivot
                   ALTERNATIVELY
                   Transpose U
                   do the same as for L
                 */
                for (j = 0; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        // check if Ui[j] is right to pivot
                        for (k = 0; k < n + 1; k++)
                        {
                            if (Up[k] <= j && j < Up[k + 1])
                            {
                                if (k > pivot)
                                {
                                    // found closest off-diagonal
                                    u_closest = k;
                                    // don't hate me
                                    goto minimum_btf;
                                    // alternatively: j=unz+1
                                }
                                break;
                            }
                        }
                    }
                }
            minimum_btf:
                pivot = MIN(l_closest, u_closest);
                // check if pivot is either already in path or n+1 (no more off-diag values)
                if (Numeric->path[pivot] == 1 || pivot == k2) // n+1)
                {
                    break;
                }
                Numeric->path[pivot] = 1;
            }
        }
    }
    // printf("Number of blocks: %d\n", nb);
    // printf("Path: ");
    // for (int i = 0; i < n - 1; i++)
    // {
    //     printf("%d, ", Numeric->path[i]);
    // }
    // printf("%d\n", Numeric->path[n-1]);
    free(Lp);
    free(Li);
    free(Lx);
    free(Up);
    free(Ui);
    free(Ux);
    free(Fi);
    free(Fp);
    free(Fx);
    free(P);
    free(Q);
    free(Qi);
    free(R);
    free(Rs);
    free(changeVector_permuted);
    return (TRUE);
}

/*
 * Secondary function. Expects a factorized matrix and "changeVector", which contains columns of A that change
 * e.g. changeVector = {3, 5}, if columns 3 and 5 contain varying entries. Needs to be permuted, since
 * L*U = P*A*Q
 */
int KLU_compute_path2(KLU_symbolic *Symbolic, KLU_numeric *Numeric, KLU_common *Common, Int *changeVector,
                      Int changeLen)
{
    /* Declarations */
    /* LU data */
    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int nb = Symbolic->nblocks;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double *Rs;

    int RET;

    /* indices and temporary variables */
    int i, k, j, ent;
    int pivot;
    int col, nextcol;
    int u_closest, l_closest;
    int flag = 1;

    // blocks
    Int k2, k1, nk;

    #ifdef FPATH
        int pivot_l, pivot_u;
    #endif

    // TODO: save sizeof(...) statically and not call for each alloc
    Int *Qi = calloc(n, sizeof(Int));
    Int *changeVector_permuted = calloc(changeLen, sizeof(Int));

    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->bpath)
    {
        KLU_free(Numeric->bpath, Numeric->nblocks, sizeof(int), Common);
    }

    Numeric->path = KLU_malloc(n, sizeof(int), Common);
    Numeric->bpath = KLU_malloc(nb, sizeof(int), Common);
    int *workpath = (int *)calloc(n, sizeof(int));

    for (i = 0; i < n; i++)
    {
        Numeric->path[i] = 0;
        workpath[i] = 0;
    }
    for (i = 0; i < nb; i++)
    {
        Numeric->bpath[i] = 0;
    }

    // TODO: solve smarter, no more klu_extracts.
    Lp = calloc(n + 1, sizeof(int));
    Up = calloc(n + 1, sizeof(int));
    Fp = calloc(n + 1, sizeof(int));
    Lx = calloc(lnz, sizeof(double));
    Ux = calloc(unz, sizeof(double));
    Fx = calloc(nzoff, sizeof(double));
    Li = calloc(lnz, sizeof(int));
    Ui = calloc(unz, sizeof(int));
    Fi = calloc(nzoff, sizeof(int));
    P = calloc(n, sizeof(int));
    Q = calloc(n, sizeof(int));
    Rs = calloc(n, sizeof(double));
    R = calloc(nb + 1, sizeof(int));

    // first, get LU decomposition
    // sloppy implementation, as there might be a smarter way to do this
    RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);

    /* check if extraction of LU matrix broke */
    if (RET != (TRUE))
    {
        return (FALSE);
    }

    // second, invert permutation vector
    // Q gives "oldcol", we need "newcol"
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    // third, apply permutation on changeVector
    for (i = 0; i < changeLen; i++)
    {
        changeVector_permuted[i] = Qi[changeVector[i]];
    }
    int* cV = (int*) calloc(n, sizeof(int));
    for ( i = 0 ; i < changeLen ; i++)
    {
        cV[changeVector_permuted[i]] = 1;
    }
    for ( i = 0; i< changeLen ; i++)
    {
        changeVector_permuted[i] = 0;
    } 
    changeLen = 0;
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            changeLen++;
        }
    }
    int ctr = 0;
    changeVector_permuted = (int*) realloc(changeVector_permuted, sizeof(int)*changeLen);
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            changeVector_permuted[ctr] = i;
            ctr++;
        }
    }
    free(cV);

    /* fourth, compute factorization path */
    if (Common->btf == FALSE)
    {
        Numeric->bpath[0] = 1;
        /* no blocks */
        for (i = 0; i < changeLen; i++)
        {
            flag = 1;
            pivot = changeVector_permuted[i];
            if (Numeric->path[pivot] == 1)
            {
                /* already computed pivot? */
                continue;
            }
            Numeric->path[pivot] = 1;
            
            #ifndef FPATH
                /* GP-based version */
                while (pivot < n && flag)
                {
                    flag = 0;
                    nextcol = Up[pivot + 1];

                    /* scan Ui for values in row "pivot" */
                    for (j = nextcol; j < unz; j++)
                    {
                        if (Ui[j] == pivot)
                        {
                            /* column j is affected
                             * find column in which j is
                            */                             
                            for (k = 0; k < n; k++)
                            {
                                if (j >= Up[k] && j < Up[k + 1])
                                {
                                    col = k;
                                    break;
                                }
                            }
                            /* col is always well-defined */
                            Numeric->path[col] = 1;
                            workpath[col] = 1;
                        }
                    }
                    for (j = pivot + 1; j < n; j++)
                    {
                        if (workpath[j] == 1)
                        {
                            workpath[j] = 0;
                            pivot = j;
                            flag = 1;
                            break;
                        }
                    }
                }
            #else
                /* off-diagonal-looking */
                pivot_l = n;
                pivot_u = n;
                while (pivot < n)
                {
                    nextcol = Lp[pivot+1];
                    /* scan Lp for next entry
                     * if this check fails, current column only has one entry
                    */
                    if(nextcol - Lp[pivot] > 1)
                    {
                        pivot_l = Li[Lp[pivot]+1];
                        /* check if found pivot is next off-diagonal entry
                         * then you can skip the next loop for U */
                        if(pivot_l - pivot == 1)
                        {
                            goto minimum_nobtf;
                        } 
                        for(j = Lp[pivot]+2; j<nextcol; j++)
                        {
                            if(Li[j]<pivot_l)
                            {
                                pivot_l = Li[j];
                            }
                        }
                    }

                    nextcol = Up[pivot + 1];
                    /* scan Ui for values in row "pivot" */
                    for (j = nextcol; j < unz; j++)
                    {
                        if (Ui[j] == pivot)
                        {
                            /* column j is affected
                             * find column in which j is
                             */
                            for (k = 0; k < n; k++)
                            {
                                if (j >= Up[k] && j < Up[k + 1])
                                {
                                    pivot_u = k;
                                    goto minimum_nobtf;
                                }
                            }
                        }
                    }
                    minimum_nobtf:
                    pivot = MIN(pivot_u, pivot_l);
                    pivot_l = n;
                    pivot_u = n;
                    Numeric->path[pivot] = 1;
                }
            #endif
        }
    }
    else
    {
        for (i = 0; i < changeLen; i++)
        {
            flag = 1;
            // get next changing column
            pivot = changeVector_permuted[i];

            // check if it was already computed
            if (Numeric->path[pivot] == 1)
            {
                // already computed pivot
                // do nothing, go to next pivot
                continue;
            }

            // set first value of singleton path
            Numeric->path[pivot] = 1;

            // find block of pivot
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */
                    Numeric->bpath[k] = 1;
                    break;
                }
            }

            if (nk == 1)
            {
                // 1x1-block, its pivot already in path
                continue;
            }

            
            #ifndef FPATH
            /* propagate until end
             * in blocks, pivot < n_block[k]
             *
             * k2 is the first column of the NEXT block
             *
             * if pivot == k2 - 1, it is the last column of the
             * current block. then, "nextcol" would be k2, which is already
             * in the next block. Thus, only do loop if pivot < k2 - 1
             */
            while (pivot < k2 - 1 && flag == 1)
            {
                flag = 0;
                nextcol = Up[pivot + 1];
                // find closest off-diagonal entries in U
                /* u saved in csc:
                   Ux = [x, x, x, x, x] entries saved column-wise
                   Ui = [0, 2, 1, 4, 3] row position of those entries in columns
                   Up = [0, 1, 2, 3, 4, 5] column pointers
                   to find closest off-diagonal of row k:
                   find all Ui[j] == pivot
                   find Up[k] <= j < Up[k+1], k > pivot
                   ALTERNATIVELY
                   Transpose U
                   do the same as for L
                 */

                /* This loop finds all nnz in row "pivot",
                 * meaning their respective columns depend on pivot column
                 */
                for (j = nextcol; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        // find column in which j is
                        for (k = k1; k < k2; k++)
                        {
                            if (j >= Up[k] && j < Up[k + 1])
                            {
                                col = k;
                                Numeric->path[col] = 1;
                                workpath[col] = 1;
                                break;
                            }
                        }
                    }
                }
                for (j = pivot + 1; j < k2; j++)
                {
                    if (workpath[j] == 1)
                    {
                        workpath[j] = 0;
                        pivot = j;
                        flag = 1;
                        break;
                    }
                }
            }
            #else
            /* propagate until end
             * in blocks, pivot < n_block[k]
             *
             * k2 is the first column of the NEXT block
             *
             * if pivot == k2 - 1, it is the last column of the
             * current block. then, "nextcol" would be k2, which is already
             * in the next block. Thus, only do loop if pivot < k2 - 1
             */
            while (pivot < k2 - 1)
            {
                pivot_l = k2;
                pivot_u = k2;

                nextcol = Lp[pivot+1];
                if(nextcol - Lp[pivot] > 1)
                {
                    pivot_l = Li[Lp[pivot]+1];
                    if(pivot_l - pivot == 1)
                    {
                        goto minimum_dobtf;
                    }
                    for(j = Lp[pivot]+2; j<nextcol; j++)
                    {
                        if(Li[j]<pivot_l)
                        {
                            pivot_l = Li[j];
                        }
                    }
                }

                nextcol = Up[pivot + 1];
                /* find closest off-diagonal entries in U */
                /* u saved in csc:
                   Ux = [x, x, x, x, x] entries saved column-wise
                   Ui = [0, 2, 1, 4, 3] row position of those entries in columns
                   Up = [0, 1, 2, 3, 4, 5] column pointers
                   to find closest off-diagonal of row k:
                   find all Ui[j] == pivot
                   find Up[k] <= j < Up[k+1], k > pivot
                   ALTERNATIVELY
                   Transpose U
                   do the same as for L
                 */

                /* This loop finds all nnz in row "pivot",
                 * meaning their respective columns depend on pivot column
                 */
                for (j = nextcol; j < unz; j++)
                {
                    if (Ui[j] == pivot)
                    {
                        // find column in which j is
                        for (k = k1; k < k2; k++)
                        {
                            if (j >= Up[k] && j < Up[k + 1])
                            {
                                pivot_u = k;
                                goto minimum_dobtf;
                            }
                        }
                    }
                }
                minimum_dobtf:
                pivot = MIN(pivot_l, pivot_u);
                Numeric->path[pivot] = 1;
                pivot_u = k2;
                pivot_l = k2;
            }
            #endif
        }
    }

    if(Common->dump == 1)
    {
        dumpAll(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, P, Q, Numeric->path, Numeric->bpath, lnz, unz, n, nzoff, nb);
    }

    free(workpath);
    free(Lp);
    free(Li);
    free(Lx);
    free(Up);
    free(Ui);
    free(Ux);
    free(Fi);
    free(Fp);
    free(Fx);
    free(P);
    free(Q);
    free(Qi);
    free(R);
    free(Rs);
    free(changeVector_permuted);
    return (TRUE);
}
