/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/* 
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"

/* 
 * "Wrapper" function. If pivot is already in path, return 1, else return 0 (i.e. pivot must be computed)
 */
int findVal(list* l, int pivot){
    int ret = search(l, pivot);
    if (ret == FOUND){
        return 1;
    }
    return 0;
}

/*
 * Main function. Expects a factorized matrix and "changeVector", which contains columns of A that change
 * e.g. changeVector = {3, 5}, if columns 3 and 5 contain varying entries. Needs to be permuted, since
 * L*U = P*A*Q
 */
int KLU_compute_path(KLU_symbolic* Symbolic, KLU_numeric* Numeric, KLU_common* Common, Int* changeVector, Int changeLen){
    // first, get LU decomposition
    // sloppy implementation, as there might be a smarter way to do this
    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double* Rs;
    int nb = Symbolic->nblocks;
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
    
    Int RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);
    if(RET != (TRUE))
        return (FALSE);

    // first and a halfly, invert permutation vector
    int* Qi = calloc(n, sizeof(int));
    for(int i=0; i<n; i++){
        Qi[Q[i]] = i;
    }

    // second, apply permutation on changeVector
    int* changeVector_permuted = calloc(sizeof(Int), changeLen);
    for(int i=0; i<changeLen; i++){
        changeVector_permuted[i] = Qi[changeVector[i]];
    }

    // third, sort permuted vector
    // in "full partial refactorisation", only find minimum value
    // sloppy selectionsort implementation for first design
    int pivot = 0;
    int indx = 0;
    for(int i=0; i<changeLen; i++){
        indx = i;
        pivot = changeVector_permuted[i];
        for(int j=i+1; j<changeLen; j++){
            if(changeVector_permuted[j]<pivot){
                pivot = changeVector_permuted[j];
                indx = j;
            }
        }
        if(indx != i){ // found smaller
            int tmp = changeVector_permuted[i];
            changeVector_permuted[i] = pivot;
            changeVector_permuted[indx] = tmp;
        }
    }

    // step two and three can / should be done externally

    // fourth, compute factorization path
    if(Numeric->path){
        free(Numeric->path);
    }
    Numeric->path = (list*)malloc(sizeof(list));
    Numeric->path->head = NULL;
    Numeric->path->tail = NULL;
    Numeric->path->length = 0;

    pivot = 0;
    int ret;
    int col, nextcol;
    int u_closest, l_closest;

    // blocks
    Int k2, k1, nk;

    for(int i=0; i<changeLen; i++){
        // get next changed column
        pivot = changeVector_permuted[i];
        
        // check if it was already computed
        ret = findVal(Numeric->path, pivot);
        if(ret == 1){
            // already computed pivot
            continue;
        }

        // singleton path, which is put at front of factorization path
        list* singleton = (list*)malloc(sizeof(list));
        singleton->head = NULL;
        singleton->tail = NULL;
        singleton->length = 0;
        
        // set first value of singleton path
        push_back(singleton, pivot);

        // find block of pivot
        for(int k=0; k<nb+1; k++){
            if(R[k] <= pivot && pivot < R[k+1]){
                k1 = R[k];
                k2 = R[k+1];
                nk = k2-k1;
                break;
            }
        }

        if(nk == 1){
            // singleton case, do nothing. pivot already added
            continue;
        }

        // propagate until end
        // in blocks, pivot < n_block[k]
        while(pivot < k2){
            u_closest = n+1;
            l_closest = n+1;

            // find closest off-diagonal entry in L
            col = Lp[pivot];
            nextcol = Lp[pivot+1];

            if(nextcol-col == 1){
                // only one entry in column => diagonal entry
                l_closest = n+1;
            } else {
                // indices are not sorted!!!! TODO. Use klu_sort maybe?

                // find closest off-diagonal entry in L, "look down" in pivot-th column
                for(int entry=col+1; entry<nextcol; entry++){
                    if(l_closest > Li[entry]){
                        l_closest = Li[entry];
                    }
                }
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
            for(int j = 0; j<unz; j++){
                if(Ui[j] == pivot){
                    // check if Ui[j] is right to pivot
                    for(int k=0; k<n+1; k++){
                        if(Up[k] <= j && j < Up[k+1]){
                            if(k > pivot){
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
            if(findVal(Numeric->path, pivot)==1||pivot==n+1)
                break;
            push_back(singleton, pivot);
        }
        push_front_list(singleton, Numeric->path);
        free(singleton);
    }
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