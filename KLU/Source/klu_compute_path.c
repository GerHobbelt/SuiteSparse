/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/* 
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"

void saveLU(double* Lx, int* Li, int* Lp, double* Ux, int* Ui, int* Up, double* Fx, int* Fi, int* Fp, int lnz, int unz, int n, int nzoff)
{
    FILE* l, *u, *g;
    g = fopen("F.txt", "w");
    l = fopen("L.txt", "w");
    u = fopen("U.txt", "w");
    int i;

    for(i=0; i<nzoff; i++)
    {
        fprintf(g, "%lf, ", Fx[i]);
    }
    fprintf(g, "\n");
    for(i=0; i<nzoff; i++)
    {
        fprintf(g, "%d, ", Fi[i]);
    }
    fprintf(g, "\n");
    for(i=0; i<n+1; i++)
    {
        fprintf(g, "%d, ", Fp[i]);
    }

    for(i=0; i<lnz; i++)
    {
        fprintf(l, "%lf, ", Lx[i]);
    }
    fprintf(l, "\n");
    for(i=0; i<lnz; i++)
    {
        fprintf(l, "%d, ", Li[i]);
    }
    fprintf(l, "\n");
    for(i=0; i<n+1; i++)
    {
        fprintf(l, "%d, ", Lp[i]);
    }

    for(i=0; i<unz; i++)
    {
        fprintf(u, "%lf, ", Ux[i]);
    }
    fprintf(u, "\n");
    for(i=0; i<unz; i++)
    {
        fprintf(u, "%d, ", Ui[i]);
    }
    fprintf(u, "\n");
    for(i=0; i<n+1; i++)
    {
        fprintf(u, "%d, ", Up[i]);
    }

    fclose(l);
    fclose(u);
    fclose(g);
}

/*
 * Main function. Expects a factorized matrix and "changeVector", which contains columns of A that change
 * e.g. changeVector = {3, 5}, if columns 3 and 5 contain varying entries. Needs to be permuted, since
 * L*U = P*A*Q
 */
int KLU_compute_path(KLU_symbolic* Symbolic, KLU_numeric* Numeric, KLU_common* Common, Int* changeVector, Int changeLen){
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
    double* Rs;
    int* Pinv = Numeric->Pinv;

    int RET;

    /* indices and temporary variables */
    int i, k, j, ent;
    int pivot;
    int col, nextcol;
    int u_closest, l_closest;

    // blocks
    Int k2, k1, nk;
    // TODO: save sizeof(...) statically and not call for each alloc
    Int* Qi = calloc(n, sizeof(Int));
    Int* changeVector_permuted = calloc(changeLen, sizeof(Int));

    Numeric->path = KLU_malloc(n, sizeof(int), Common);
    Numeric->bpath = KLU_malloc(nb, sizeof(int), Common);

    for(i = 0 ; i < n ; i++)
    {
        Numeric->path[i] = 0;
    }
    for(i = 0 ; i < nb ; i++)
    {
        Numeric->bpath[i] = 0;
    }

    // TODO: solve smarter, no more klu_extracts.
    Lp = calloc(n+1, sizeof(int));
    Up = calloc(n+1, sizeof(int));
    Fp = calloc(n+1, sizeof(int));
    Lx = calloc(lnz, sizeof(double));
    Ux = calloc(unz, sizeof(double));
    Fx = calloc(nzoff, sizeof(double));
    Li = calloc(lnz, sizeof(int));
    Ui = calloc(unz, sizeof(int));
    Fi = calloc(nzoff, sizeof(int));
    P = calloc(n, sizeof(int));
    Q = calloc(n, sizeof(int));
    Rs = calloc(n, sizeof(double));
    R = calloc(nb+1, sizeof(int));
    
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
    for (i = 0 ; i < n ; i++)
    {
        Qi [Q [i]] = i;
    }

    // third, apply permutation on changeVector
    for (i = 0 ; i < changeLen ; i++)
    {
        changeVector_permuted [i] = Qi [changeVector [i]];
    }

    // fourth, sort permuted vector
    // in "full partial refactorisation", only find minimum value
    // sloppy selectionsort implementation for first design
    for (i = 0 ; i < changeLen ; i++)
    {
        k = i;
        pivot = changeVector_permuted [i];
        for (j = i + 1 ; j < changeLen ; j++)
        {
            if (changeVector_permuted [j] < pivot)
            {
                pivot = changeVector_permuted [j];
                k = j;
            }
        }
        if (k != i)
        { 
            // found smaller
            // switch positions
            int tmp = changeVector_permuted [i];
            changeVector_permuted [i] = pivot;
            changeVector_permuted [k] = tmp;
        }
    }

    // step three and four can / should be done externally

    // fifth, compute factorization path
    for (i = 0 ; i < changeLen ; i++)
    {
        // get next changing column
        pivot = changeVector_permuted [i];
        
        // check if it was already computed
        if (Numeric->path[pivot] == 1)
        {
            // already computed pivot
            // do nothing, go to next pivot
            continue;
        }

        // singleton path, which is put at front of factorization path
        // list* singleton = (list*)malloc(sizeof(list));
        // singleton->head = NULL;
        // singleton->tail = NULL;
        // singleton->length = 0;
        
        // set first value of singleton path
        // push_back(singleton, pivot);
        Numeric->path[pivot] = 1;

        // find block of pivot
        for (k = 0 ; k < nb+1 ; k++)
        {
            if (R [k] <= pivot && pivot < R [k+1])
            {
                k1 = R [k];
                /* set varying block */
                Numeric->bpath[k] = 1;
                k2 = R [k+1];
                nk = k2 - k1;
                break;
            }
        }

        if (nk == 1)
        {
            // Maybe redundant?
            // singleton case, put pivot in front of path
            // push_front_list(singleton, Numeric->path);
            // free(singleton);
            continue;
        }

        // propagate until end
        // in blocks, pivot < n_block[k]
        while (pivot < k2)
        {
            //u_closest = n+1;
            //l_closest = n+1;
            u_closest = k2;
            l_closest = k2;

            // find closest off-diagonal entry in L
            col = Lp[pivot];
            nextcol = Lp[pivot+1];

            if (nextcol - col == 1)
            {
                // only one entry in column => diagonal entry
                //l_closest = n+1;
                l_closest = k2;
            } 
            else 
            {
                // indices are not sorted!!!! TODO. Use klu_sort maybe?

                // find closest off-diagonal entry in L, "look down" in pivot-th column
                for (ent = col + 1 ; ent < nextcol ; ent++)
                {
                    if (l_closest > Li [ent])
                    {
                        l_closest = Li [ent];
                    }
                }
            }
            if(l_closest - pivot == 1)
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
            for (j = 0 ; j < unz ; j++)
            {
                if(Ui [j] == pivot)
                {
                    // check if Ui[j] is right to pivot
                    for (k = 0 ; k < n+1 ; k++)
                    {
                        if (Up [k] <= j && j < Up [k+1])
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
            if (Numeric->path[pivot] == 1 || pivot == k2)//n+1)
            {
                break;
            }
            //push_back(singleton, pivot);
            Numeric->path[pivot] = 1;
        }
        // push_front_list(singleton, Numeric->path);
        // free(singleton);
    }
    //sort_list(Numeric->path, QUICK);
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