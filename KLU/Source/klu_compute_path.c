/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/*
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"

int KLU_extract_quick(
    /* inputs: */
    KLU_numeric *Numeric,
    KLU_symbolic *Symbolic,

    /* outputs: */
    Int *Ui,
    Int *Up,
    Int *Q,
    Int *R
)
{
    /* placeholder */
    Entry *Lx2, *Ux2, *Ukk ;
    Unit* LU;
    Int *Lip, *Llen, *Uip, *Ulen, *Li2, *Ui2 ;
    Int nz, k1, k2, block, nk, kk, len, p, n, nblocks, k;
    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;

    /* ---------------------------------------------------------------------- */
    /* extract block boundaries */
    /* ---------------------------------------------------------------------- */

    if (R != NULL)
    {
        for (block = 0 ; block <= nblocks ; block++)
        {
            R [block] = Symbolic->R [block] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* extract column permutation */
    /* ---------------------------------------------------------------------- */

    if (Q != NULL)
    {
        for (k = 0 ; k < n ; k++)
        {
            Q [k] = Symbolic->Q [k] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* extract each block of U */
    /* ---------------------------------------------------------------------- */

    if (Up != NULL && Ui != NULL)
    {
        nz = 0 ;
        for (block = 0 ; block < nblocks ; block++)
        {
            k1 = Symbolic->R [block] ;
            k2 = Symbolic->R [block+1] ;
            nk = k2 - k1 ;
            if (nk == 1)
            {
                /* singleton block */
                Up [k1] = nz ;
                Ui [nz] = k1 ;
                nz++ ;
            }
            else
            {
                /* non-singleton block */
                LU = Numeric->LUbx [block] ;
                Uip = Numeric->Uip + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                for (kk = 0 ; kk < nk ; kk++)
                {
                    Up [k1+kk] = nz ;
                    GET_POINTER (LU, Uip, Ulen, Ui2, Ux2, kk, len) ;
                    for (p = 0 ; p < len ; p++)
                    {
                        Ui [nz] = k1 + Ui2 [p] ;
                        nz++ ;
                    }
                    /* add the diagonal entry */
                    Ui [nz] = k1 + kk ;
                    nz++ ;
                }
            }
        }
        Up [n] = nz ;
        ASSERT (nz == Numeric->unz) ;
    }
    return (TRUE);
}

/*
 * Secondary function. Expects a factorized matrix and "changeVector", which contains columns of A that change
 * e.g. changeVector = {3, 5}, if columns 3 and 5 contain varying entries. Needs to be permuted, since
 * L*U = P*A*Q
 */
int KLU_compute_path(
                    KLU_symbolic *Symbolic, 
                    KLU_numeric *Numeric, 
                    KLU_common *Common, 
                    Int *changeVector,
                    Int changeLen
                    )
{
    /* This function is very long, because you have to implement BTF and no BTF-case... */
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
    }
    for (i = 0; i < nb; i++)
    {
        Numeric->bpath[i] = 0;
    }


    // Up = calloc(n + 1, sizeof(int));
    // Ui = calloc(unz, sizeof(int));
    // Q = calloc(n, sizeof(int));
    // R = calloc(nb+1, sizeof(int));

    // // first, get LU decomposition
    // // sloppy implementation, as there might be a smarter way to do this
    // RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);
    //RET = KLU_extract_quick(Numeric, Symbolic, Ui, Up, Q, R);
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
        }
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
