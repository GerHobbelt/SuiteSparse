/* ========================================================================== */
/* === KLU_compute_path ===================================================== */
/* ========================================================================== */

/*
 * Computes the factorization path given change vector
 * Any new refactorisation can be computed by iterating over entries in factorization path
 * instead of all columns (see klu_partial.c)
 */
#include "klu_internal.h"

// int KLU_extract_quick(
//     /* inputs: */
//     KLU_numeric *Numeric,
//     KLU_symbolic *Symbolic,

//     /* outputs: */
//     Int *Ui,
//     Int *Up,
//     Int *Q,
//     Int *R
// )
// {
//     /* placeholder */
//     Entry *Lx2, *Ux2, *Ukk ;
//     Unit* LU;
//     Int *Lip, *Llen, *Uip, *Ulen, *Li2, *Ui2 ;
//     Int nz, k1, k2, block, nk, kk, len, p, n, nblocks, k;
//     n = Symbolic->n ;
//     nblocks = Symbolic->nblocks ;

//     /* ---------------------------------------------------------------------- */
//     /* extract block boundaries */
//     /* ---------------------------------------------------------------------- */

//     if (R != NULL)
//     {
//         for (block = 0 ; block <= nblocks ; block++)
//         {
//             R [block] = Symbolic->R [block] ;
//         }
//     }

//     /* ---------------------------------------------------------------------- */
//     /* extract column permutation */
//     /* ---------------------------------------------------------------------- */

//     if (Q != NULL)
//     {
//         for (k = 0 ; k < n ; k++)
//         {
//             Q [k] = Symbolic->Q [k] ;
//         }
//     }

//     /* ---------------------------------------------------------------------- */
//     /* extract each block of U */
//     /* ---------------------------------------------------------------------- */

//     if (Up != NULL && Ui != NULL)
//     {
//         nz = 0 ;
//         for (block = 0 ; block < nblocks ; block++)
//         {
//             k1 = Symbolic->R [block] ;
//             k2 = Symbolic->R [block+1] ;
//             nk = k2 - k1 ;
//             if (nk == 1)
//             {
//                 /* singleton block */
//                 Up [k1] = nz ;
//                 Ui [nz] = k1 ;
//                 nz++ ;
//             }
//             else
//             {
//                 /* non-singleton block */
//                 LU = Numeric->LUbx [block] ;
//                 Uip = Numeric->Uip + k1 ;
//                 Ulen = Numeric->Ulen + k1 ;
//                 for (kk = 0 ; kk < nk ; kk++)
//                 {
//                     Up [k1+kk] = nz ;
//                     GET_POINTER (LU, Uip, Ulen, Ui2, Ux2, kk, len) ;
//                     for (p = 0 ; p < len ; p++)
//                     {
//                         Ui [nz] = k1 + Ui2 [p] ;
//                         nz++ ;
//                     }
//                     /* add the diagonal entry */
//                     Ui [nz] = k1 + kk ;
//                     nz++ ;
//                 }
//             }
//         }
//         Up [n] = nz ;
//         ASSERT (nz == Numeric->unz) ;
//     }
//     return (TRUE);
// }

/*
 * Computes Factorization Path.
 */
int KLU_compute_path(
                    KLU_symbolic *Symbolic, 
                    KLU_numeric *Numeric, 
                    KLU_common *Common, 
                    Int Ap [ ],
                    Int Ai [ ],
                    Int *variable_columns,
                    Int *variable_rows,
                    Int n_variable_entries
                    )
{

    /* This method computes the factorization path.
     * "output": 
     *          - variable_offdiag_orig_entry: position of entries in Ax, which end up in off-diagonal block F
     *          - variable_offdiag_perm_entry: position of entries in F, which are varying in Ax
     *          - bpath: array of length nblocks+1. bpath[k]...bpath[k+1]-1 indicate the variable columns in block k
     *          - path: path[bpath[k]]...path[bpath[k+1]-1] contain the variable columns in block k
     *      
     */

    if(variable_columns == NULL)
    {
        return FALSE;
    }
    if(variable_rows == NULL)
    {
        return FALSE;
    }
    if(n_variable_entries <= 0)
    {
        return TRUE;
    }
    /* This function is very long, because you have to implement BTF and no BTF-case... */
    /* Declarations */
    /* LU data */

    int n = Symbolic->n;
    int lnz = Numeric->lnz;
    int unz = Numeric->unz;
    int nzoff = Numeric->nzoff;
    int nb = Symbolic->nblocks;
    int *Pinv = Numeric->Pinv;
    int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
    double *Lx, *Ux, *Fx;
    int *P, *Q, *R;
    double *Rs;
    int RET;
    int oldcol, pend, p, newrow;
    int poff = 0, ctr = 0, variable_offdiag_length = 0;

    Int n_variable_entries_new = n_variable_entries;

    Int* variable_offdiag_orig_entry = (Int*)calloc(nzoff, sizeof(Int));
    Int *variable_offdiag_perm_entry = (Int*)calloc(nzoff, sizeof(Int));

    /* indices and temporary variables */
    int i, k, j, ent, block;
    int pivot;
    int col, nextcol;
    int flag = 1;

    /* blocks */
    Int k2, k1, nk;

    Int *Qi = calloc(n, sizeof(Int));
    Int *variable_columns_in_LU = calloc(n_variable_entries, sizeof(Int));
    Int *variable_rows_in_LU = calloc(n_variable_entries, sizeof(Int));

    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->block_path)
    {
        KLU_free(Numeric->block_path, Numeric->nblocks, sizeof(int), Common);
    }
    if (Numeric->variable_block)
    {
        KLU_free(Numeric->variable_block, Numeric->n_variable_blocks, sizeof(int), Common);
    }
    if (Numeric->variable_offdiag_orig_entry)
    {
        KLU_free(Numeric->variable_offdiag_orig_entry, Numeric->variable_offdiag_length, sizeof(int), Common);
    }
    if (Numeric->variable_offdiag_perm_entry)
    {
        KLU_free(Numeric->variable_offdiag_perm_entry, Numeric->variable_offdiag_length, sizeof(int), Common);
    }

    Numeric->block_path = KLU_malloc(nb, sizeof(int), Common);
    int workpath[n];
    Int path[n];
    Int nvblocks[nb];
    
    for (i  = 0 ; i < n ; i++)
    {
        path[i] = 0;
    }

    for (i = 0; i < nb; i++)
    {
        nvblocks[i] = 0;
        Numeric->block_path[i] = 0;
    }

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

    /* ---------------------------------------------------------------- */
    /* first, get LU decomposition */
    /* sloppy implementation, as there might be a smarter way to do this */
    /* ---------------------------------------------------------------- */
    RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);

    /* check if extraction of LU matrix broke */
    if (RET != (TRUE))
    {
        return (FALSE);
    }

    /* ---------------------------------------------------------------- */
    /* second, invert permutation vector */
    /* Q gives "oldcol", we need "newcol" */
    /* ---------------------------------------------------------------- */
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    /* ---------------------------------------------------------------- */
    /* third, apply permutation on variable_columns */
    /* ---------------------------------------------------------------- */

    for (i = 0; i < n_variable_entries; i++)
    {
        variable_columns_in_LU[i] = Qi[variable_columns[i]];
        variable_rows_in_LU[i] = Pinv[variable_rows[i]];
    }

    /* ---------------------------------------------------------------- */
    /* fourth, determine variable off-diagonal entries */
    /* only necessary if BTF used */
    /* if BTF is not used, no off-diagonal blocks are available */
    /* ---------------------------------------------------------------- */

    if(Common->btf == TRUE)
    {
        /* iterate over all blocks */
        for(block = 0 ; block < nb ; block++)
        {
            k1 = R [block];
            k2 = R [block+1];
            nk = k2 - k1;
            for (k = 0 ; k < nk ; k++)
            {
                oldcol = Q [k+k1] ;
                pend = Ap [oldcol+1] ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    newrow = Pinv [Ai [p]] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */
                        poff++ ;

                        /* check if entry is in variable column */
                        /* TODO */
                        for(i = 0 ; i < n_variable_entries ; i++)
                        {
                            if(variable_columns_in_LU[i] == k+k1 && variable_rows_in_LU[i] == newrow)
                            {
                                /* set to -1, because they're off-diagonal now */
                                variable_columns_in_LU[i] = -1;
                                variable_rows_in_LU[i] = -1;

                                variable_offdiag_orig_entry[variable_offdiag_length] = p;
                                variable_offdiag_perm_entry[variable_offdiag_length++] = poff;

                                n_variable_entries_new--;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    ASSERT(variable_offdiag_length == n_variable_entries - n_variable_entries_new);

    Numeric->variable_offdiag_orig_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_perm_entry = KLU_malloc(variable_offdiag_length, sizeof(Int), Common);
    Numeric->variable_offdiag_length = variable_offdiag_length;

    for(i = 0; i < variable_offdiag_length ; i++)
    {
        Numeric->variable_offdiag_orig_entry[i] = variable_offdiag_orig_entry[i];
        Numeric->variable_offdiag_perm_entry[i] = variable_offdiag_perm_entry[i];
    }

    int* cV = calloc(n, sizeof(int));
    for(i = 0; i < n_variable_entries ; i++)
    {
        if(variable_columns_in_LU[i] != -1)
        {
            cV[variable_columns_in_LU[i]] = 1;
        }
    }
    for( i = 0; i < n_variable_entries ; i++)
    {
        variable_columns_in_LU[i] = 0;
    }
    variable_columns_in_LU = (Int*)realloc(variable_columns_in_LU, sizeof(Int)*n_variable_entries_new);
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] != 0)
        {
            variable_columns_in_LU[ctr++] = i;
        }
    }
    free(cV);

    ASSERT(ctr == variable_entries_new);

    /* ---------------------------------------------------------------- */
    /* fifth, compute factorization path of blocks */
    /* ---------------------------------------------------------------- */

    if (Common->btf == FALSE)
    {
        Numeric->variable_block = KLU_malloc(1, sizeof(int), Common);
        Numeric->block_path[0] = 0;
        Numeric->variable_block[0] = 0;
        Numeric->n_variable_blocks = 1;
        /* no blocks */
        for (i = 0; i < n_variable_entries_new; i++)
        {
            flag = 1;
            pivot = variable_columns_in_LU[i];
            if (path[pivot] == 1)
            {
                /* already computed pivot? */
                continue;
            }
            path[pivot] = 1;
            
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
                        path[col] = 1;
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
        ctr = 0;
        for( i = 0 ; i < n ; i++)
        {
            if(path[i] == 1)
            {
                ctr++;
            }
        }
        Numeric->path = KLU_malloc(ctr, sizeof(int), Common);
        for(i = 0; i < ctr ; i++)
        {
            Numeric->path[i] = 0;
        }
        Numeric->block_path[1] = ctr+1;
        ctr = 0;
        for(i = 0; i < n ; i++)
        {
            if(path[i] == 1)
            {
                Numeric->path[ctr++] = i;
            }
        }
    }
    else
    {
        Numeric->n_variable_blocks = 0;
        for (i = 0; i < n_variable_entries; i++)
        {
            flag = 1;
            /* get next changing column */
            pivot = variable_columns_in_LU[i];

            /* check if it was already computed */
            if (path[pivot] == 1)
            {
                /* already computed pivot
                 * do nothing, go to next pivot
                 */
                continue;
            }

            /* set first value of singleton path */
            path[pivot] = 1;

            /* find block of pivot */
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */

                    if(nvblocks[k] != 1)
                    {
                        nvblocks[k] = 1;
                        Numeric->n_variable_blocks += 1;
                    }
                    break;
                }
            }

            if (nk == 1)
            {
                /* 1x1-block, its pivot already in path */
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
                                col = k;
                                path[col] = 1;
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

        /* count how many variable columns there are */
        ctr = 0;
        for( i = 0 ; i < n ; i++)
        {
            if(path[i] == 1)
            {
                ctr++;
            }
        }

        /* allocate and initialize memory */
        Numeric->path = KLU_malloc(ctr, sizeof(int), Common);
        Numeric->variable_block = KLU_malloc(Numeric->n_variable_blocks, sizeof(int), Common);
        for(i = 0; i < ctr ; i++)
        {
            Numeric->path[i] = 0;
        }

        /* set variable blocks */
        ctr = 0;
        for(i = 0 ; i < nb ; i++)
        {
            if(nvblocks[i] == 1)
            {
                Numeric->variable_block[ctr++] = i;
            }
        }

        ASSERT(ctr == Numeric->n_variable_blocks);

        /* assemble factorization path */
        ctr = 0;
        for(i = 0; i < n ; i++)
        {
            if(path[i] == 1)
            {
                Numeric->path[ctr++] = i;
            }
        }
        block = 0;

        /* determine from where to where in each block there are variable columns */
        k = 0;
        for(i = 0 ; i < Numeric->n_variable_blocks ; i++)
        {
            k1 = R[Numeric->variable_block[i]]; 
            k2 = R[Numeric->variable_block[i]+1];
            Numeric->block_path[Numeric->variable_block[i]] = k;
            while(k < ctr && Numeric->path[k] < k2 - 1)
            {
                k++;
            }
            k++;
            Numeric->block_path[Numeric->variable_block[i]+1] = k;
        }
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
    free(variable_columns_in_LU);
    free(variable_rows_in_LU);
    return (TRUE);
}

/*
 * This function determines the first varying column for 
 * partial refactorization by refactorization restart (PR-RR)
 */
int KLU_determine_start(
        KLU_symbolic *Symbolic, 
        KLU_numeric *Numeric, 
        KLU_common *Common, 
        Int *variable_columns,
        Int n_variable_entries
    )
{
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
    int i, k, j, ent, block;
    int pivot;
    int col, nextcol;
    int u_closest, l_closest;
    int flag = 1;

    /* blocks */
    Int k2, k1, nk;

    /* TODO: save sizeof(...) statically and not call for each alloc */
    Int *Qi = calloc(n, sizeof(Int));
    Int *variable_columns_in_LU = calloc(n_variable_entries, sizeof(Int));

    if (Numeric->path)
    {
        KLU_free(Numeric->path, n, sizeof(int), Common);
    }
    if (Numeric->block_path)
    {
        KLU_free(Numeric->block_path, Numeric->nblocks, sizeof(int), Common);
    }
    if (Numeric->start)
    {
        KLU_free(Numeric->start, Numeric->nblocks, sizeof(int), Common);
    }

    /* Numeric->path = KLU_malloc(n, sizeof(int), Common); */
    Numeric->block_path = KLU_malloc(nb, sizeof(int), Common);
    Numeric->start = KLU_malloc(nb, sizeof(int), Common);
/*
    for (i = 0; i < n; i++)
    {
        Numeric->path[i] = 0;
    }
*/
    for (i = 0; i < nb; i++)
    {
        Numeric->block_path[i] = 0;
        Numeric->start[i] = n;
    }

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

    /* first, get LU decomposition
     * sloppy implementation, as there might be a smarter way to do this
     */
    RET = klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);

    /* check if extraction of LU matrix broke */
    if (RET != (TRUE))
    {
        return (FALSE);
    }

    /* second, invert permutation vector
     * Q gives "oldcol", we need "newcol"
     */
    for (i = 0; i < n; i++)
    {
        Qi[Q[i]] = i;
    }

    /* third, apply permutation on variable_columns */
    for (i = 0; i < n_variable_entries; i++)
    {
        variable_columns_in_LU[i] = Qi[variable_columns[i]];
    }
    int* cV = (int*) calloc(n, sizeof(int));
    for ( i = 0 ; i < n_variable_entries ; i++)
    {
        cV[variable_columns_in_LU[i]] = 1;
    }
    for ( i = 0; i< n_variable_entries ; i++)
    {
        variable_columns_in_LU[i] = 0;
    } 
    n_variable_entries = 0;
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            n_variable_entries++;
        }
    }
    int ctr = 0;
    variable_columns_in_LU = (int*) realloc(variable_columns_in_LU, sizeof(int)*n_variable_entries);
    for ( i = 0 ; i < n ; i++)
    {
        if(cV[i] == 1)
        {
            variable_columns_in_LU[ctr] = i;
            ctr++;
        }
    }
    free(cV);

    if (Common->btf == FALSE)
    {
        /* no btf case
         * => identify first varying entry in entire matrix */
        Numeric->block_path[0] = 0;
        Numeric->block_path[1] = 1;
        /* Numeric->variable_block[0] = 0; */

        /* find minimum in variable_columns_in_LU */
        pivot = variable_columns_in_LU[0];
        for(i = 1; i < n ; i++)
        {
            if(pivot > variable_columns_in_LU[i])
            {
                pivot = variable_columns_in_LU[i];
            }
        }

        /* set first varying column internally */
        Numeric->start[0] = pivot;
    }
    else
    {
        /* btf case
         * 
         * iterate over variable_columns_in_LU
         * for each variable column
         * find its block
         *      put block in block-path
         *      if block larger than 1
         *          start[block] = min(start[block], column), if start[block] != 0
         */

         for(i = 0; i < n_variable_entries ; i++)
         {
            /* grab variable column number i */
            pivot = variable_columns_in_LU[i];

            /* find block */
            for (k = 0; k < nb; k++)
            {
                k1 = R[k];
                k2 = R[k + 1];
                if (k1 <= pivot && pivot < k2)
                {
                    nk = k2 - k1;

                    /* set varying block */
                    /* TODO!!! */
                    /* Numeric->block_path[k] = 1; */
                    break;
                }
            }

            if (nk == 1)
            {
                /* 1x1-block, nothing left to do */
                continue;
            }

            /* nk x nk block. check if pivot column is smaller (before) current minimum */
            if(Numeric->start[k] > pivot)
            {
                Numeric->start[k] = pivot;
            }
         }
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
    free(variable_columns_in_LU);
    return (TRUE);
}