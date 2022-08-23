/* ========================================================================== */
/* === KLU_partial ========================================================== */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with KLU_analyze,
 * factoring it once with KLU_factor, and computing factorization path.
 * This routine cannot do any numerical pivoting.  The pattern of the
 * input matrix (Ap, Ai) must be identical to the pattern given to
 * KLU_factor.
 */

#include "klu_internal.h"

/* ========================================================================== */
/* === KLU_partial_refactorization_restart  ================================= */
/* ========================================================================== */

Int KLU_partial_refactorization_restart /* returns TRUE if successful, FALSE otherwise */
    (
        /* inputs, not modified */
        Int Ap[],                            /* size n+1, column pointers */
        Int Ai[],                            /* size nz, row indices */
        double Ax[], 
        KLU_symbolic *Symbolic, /* now also contains factorization path */
        /* input/output */
        KLU_numeric *Numeric, 
        KLU_common *Common)
{
    Entry ukk, ujk, s;
    Entry *Offx, *Lx, *Ux, *X, *Az, *Udiag;
    double *Rs;
    Int *Q, *R, *Pnum, *Ui, *Li, *Pinv, *Lip, *Uip, *Llen, *Ulen;
    Unit **LUbx;
    Unit *LU;
    Int k1, k2, nk, k, block, oldcol, pend, oldrow, n, p, newrow, scale, nblocks, poff, i, j, up, ulen, llen, maxblock,
        nzoff;

    Int z = 0;
    Int doRefact = 0;

    Int* start = Numeric->start;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE);
    }
    Common->status = KLU_OK;

    if (Numeric == NULL)
    {
        /* invalid Numeric object */
        Common->status = KLU_INVALID;
        return (FALSE);
    }

    Common->numerical_rank = EMPTY;
    Common->singular_col = EMPTY;

    Az = (Entry *)Ax;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n;
    Q = Symbolic->Q;
    R = Symbolic->R;
    nblocks = Symbolic->nblocks;
    maxblock = Symbolic->maxblock;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    Pnum = Numeric->Pnum;
    Offx = (Entry *)Numeric->Offx;

    LUbx = (Unit **)Numeric->LUbx;

    scale = Common->scale;
    if (scale > 0)
    {
        /* factorization was not scaled, but refactorization is scaled */
        if (Numeric->Rs == NULL)
        {
            Numeric->Rs = KLU_malloc(n, sizeof(double), Common);
            if (Common->status < KLU_OK)
            {
                Common->status = KLU_OUT_OF_MEMORY;
                return (FALSE);
            }
        }
    }
    else
    {
        /* no scaling for refactorization; ensure Numeric->Rs is freed.  This
         * does nothing if Numeric->Rs is already NULL. */
        Numeric->Rs = KLU_free(Numeric->Rs, n, sizeof(double), Common);
    }
    Rs = Numeric->Rs;

    Pinv = Numeric->Pinv;
    X = (Entry *)Numeric->Xwork;
    Common->nrealloc = 0;
    Udiag = Numeric->Udiag;
    nzoff = Symbolic->nzoff;
    /* ---------------------------------------------------------------------- */
    /* check the input matrix compute the row scale factors, Rs */
    /* ---------------------------------------------------------------------- */

    /* permute scaling Rs back */
    if (scale > 0)
    {
        for (k = 0; k < n; k++)
        {
            REAL(X[k]) = Rs[Pinv[k]];
        }
        for (k = 0; k < n; k++)
        {
            Rs[k] = REAL(X[k]);
        }
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace X */
    /* ---------------------------------------------------------------------- */

    for (k = 0; k < maxblock; k++)
    {
        /* X [k] = 0 */
        CLEAR(X[k]);
    }

    /* ---------------------------------------------------------------------- */
    /* assemble off-diagonal blocks */
    /* ---------------------------------------------------------------------- */

    poff = 0;

    /* ---------------------------------------------------------------------- */
    /* factor each block */
    /* ---------------------------------------------------------------------- */

    if (scale <= 0)
    {
        /* ------------------------------------------------------------------ */
        /* no scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0; block < nblocks; block++)
        {
            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R[block];
            k2 = R[block + 1];
            nk = k2 - k1;

            if (nk == 1)
            {
                /* ---------------------------------------------------------- */
                /* singleton case */
                /* ---------------------------------------------------------- */

                if (Numeric->block_path[block] != 1)
                {
                    /* ----------------------------------------------------------
                     */
                    /* no refactorization, only raise block index ctr. */
                    /* ----------------------------------------------------------
                     */

                    oldcol = Q[k1];
                    pend = Ap[oldcol + 1];
                    for (p = Ap[oldcol]; p < pend; p++)
                    {
                        newrow = Pinv[Ai[p]] - k1;
                        if (newrow < 0 && poff < nzoff)
                        {
                            poff++;
                        }
                    }
                }
                else
                {
                    oldcol = Q[k1];
                    pend = Ap[oldcol + 1];
                    CLEAR(s);
                    for (p = Ap[oldcol]; p < pend; p++)
                    {
                        newrow = Pinv[Ai[p]] - k1;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal block */
                            Offx[poff] = Az[p];
                            poff++;
                        }
                        else
                        {
                            /* singleton */
                            s = Az[p];
                        }
                    }
                    Udiag[k1] = s;
                }
            }
            else
            {
                if (Numeric->block_path[block] != 1)
                {
                    /* 
                      unlikely
                      encountered a block > 1 that has no refact. effort, i.e. no varying entries
                      only raise counter 
                    */
                    for (k = 0; k < nk; k++)
                    {
                        oldcol = Q[k + k1];
                        pend = Ap[oldcol + 1];
                        for (p = Ap[oldcol]; p < pend; p++)
                        {
                            newrow = Pinv[Ai[p]] - k1;
                            if (newrow < 0 && poff < nzoff)
                            {
                                poff++;
                            }
                        }
                    }
                }
                else
                {
                    /* ----------------------------------------------------------
                     */
                    /* construct and factor the kth block */
                    /* ----------------------------------------------------------
                     */
                    Lip = Numeric->Lip + k1;
                    Llen = Numeric->Llen + k1;
                    Uip = Numeric->Uip + k1;
                    Ulen = Numeric->Ulen + k1;
                    LU = LUbx[block];

                    for (k = 0; k < nk; k++)
                    {

                        if(k+k1 < start[block])
                        {
                            /* block contains varying entries, but column k of
                             * block has no refactorization effort only raise
                             * counter */
                            oldcol = Q[k + k1];
                            pend = Ap[oldcol + 1];
                            for (p = Ap[oldcol]; p < pend; p++)
                            {
                                newrow = Pinv[Ai[p]] - k1;
                                if (newrow < 0 && poff < nzoff)
                                {
                                    poff++;
                                }
                            }
                        }
                        else
                        {
                            /* ------------------------------------------------------
                            */
                            /* scatter kth column of the block into workspace X
                            */
                            /* ------------------------------------------------------
                            */

                            oldcol = Q[k + k1];
                            pend = Ap[oldcol + 1];
                            for (p = Ap[oldcol]; p < pend; p++)
                            {
                                newrow = Pinv[Ai[p]] - k1;
                                if (newrow < 0 && poff < nzoff)
                                {
                                    /* entry in off-diagonal block */
                                    Offx[poff] = Az[p];
                                    poff++;
                                }
                                else
                                {
                                    /* (newrow,k) is an entry in the block */
                                    X[newrow] = Az[p];
                                }
                            }

                            /* ------------------------------------------------------
                            */
                            /* compute kth column of U, and update kth column of
                            * A */
                            /* ------------------------------------------------------
                            */

                            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, ulen);
                            for (up = 0; up < ulen; up++)
                            {
                                j = Ui[up];
                                ujk = X[j];
                                /* X [j] = 0 */
                                CLEAR(X[j]);
                                Ux[up] = ujk;
                                GET_POINTER(LU, Lip, Llen, Li, Lx, j, llen);
                                for (p = 0; p < llen; p++)
                                {
                                    /* X [Li [p]] -= Lx [p] * ujk */
                                    MULT_SUB(X[Li[p]], Lx[p], ujk);
                                }
                            }
                            /* get the diagonal entry of U */
                            ukk = X[k];
                            /* X [k] = 0 */
                            CLEAR(X[k]);
                            /* singular case */
                            if (IS_ZERO(ukk))
                            {
                                /* matrix is numerically singular */
                                Common->status = KLU_SINGULAR;
                                if (Common->numerical_rank == EMPTY)
                                {
                                    Common->numerical_rank = k + k1;
                                    Common->singular_col = Q[k + k1];
                                }
                                if (Common->halt_if_singular)
                                {
                                    /* do not continue the factorization */
                                    return (FALSE);
                                }
                            }
                            /* pivot vadility testing */
                            else if (SCALAR_ABS(ukk) < Common->pivot_tol_fail)
                            {
                                /* pivot is too small */
                                Common->status = KLU_PIVOT_FAULT;
                                if (Common->halt_if_pivot_fails)
                                {
                                    /* do not continue the factorization */
                                    return (FALSE);
                                }
                            }
                            Udiag[k + k1] = ukk;
                            /* gather and divide by pivot to get kth column of L
                            */
                            GET_POINTER(LU, Lip, Llen, Li, Lx, k, llen);

                            for (p = 0; p < llen; p++)
                            {
                                i = Li[p];
                                DIV(Lx[p], X[i], ukk);
                                CLEAR(X[i]);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        /* ------------------------------------------------------------------ */
        /* scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0; block < nblocks; block++)
        {
            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R[block];
            k2 = R[block + 1];
            nk = k2 - k1;

            if (nk == 1)
            {
                /* ---------------------------------------------------------- */
                /* singleton case */
                /* ---------------------------------------------------------- */
                if (Numeric->block_path[block] != 1)
                {
                    /* --------------------------------------------------------- */
                    /* no refactorization, only raise block index ctr. */
                    /* --------------------------------------------------------- */

                    oldcol = Q[k1];
                    pend = Ap[oldcol + 1];
                    for (p = Ap[oldcol]; p < pend; p++)
                    {
                        newrow = Pinv[Ai[p]] - k1;
                        if (newrow < 0 && poff < nzoff)
                        {
                            poff++;
                        }
                    }
                }
                else
                {
                    oldcol = Q[k1];
                    pend = Ap[oldcol + 1];
                    CLEAR(s);
                    for (p = Ap[oldcol]; p < pend; p++)
                    {
                        oldrow = Ai[p];
                        newrow = Pinv[oldrow] - k1;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal block */
                            /* Offx [poff] = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN(Offx[poff], Az[p], Rs[oldrow]);
                            poff++;
                        }
                        else
                        {
                            /* singleton */
                            /* s = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN(s, Az[p], Rs[oldrow]);
                        }
                    }
                    Udiag[k1] = s;
                }
            }
            else
            {
                /* ---------------------------------------------------------- */
                /* construct and factor the kth block */
                /* ---------------------------------------------------------- */
                if (Numeric->block_path[block] != 1)
                {
                    // unlikely
                    // encountered a block > 1 that has no refact. effort
                    // only raise counter
                    for (k = 0; k < nk; k++)
                    {
                        oldcol = Q[k + k1];
                        pend = Ap[oldcol + 1];
                        for (p = Ap[oldcol]; p < pend; p++)
                        {
                            newrow = Pinv[Ai[p]] - k1;
                            if (newrow < 0 && poff < nzoff)
                            {
                                poff++;
                            }
                        }
                    }
                }
                else
                {
                    Lip = Numeric->Lip + k1;
                    Llen = Numeric->Llen + k1;
                    Uip = Numeric->Uip + k1;
                    Ulen = Numeric->Ulen + k1;
                    LU = LUbx[block];

                    for (k = 0; k < nk; k++)
                    {
                        if(k+k1 < start[block])
                        {
                            /* block contains varying entries, but column k of
                             * block has no refactorization effort only raise
                             * counter */
                            oldcol = Q[k + k1];
                            pend = Ap[oldcol + 1];
                            for (p = Ap[oldcol]; p < pend; p++)
                            {
                                newrow = Pinv[Ai[p]] - k1;
                                if (newrow < 0 && poff < nzoff)
                                {
                                    poff++;
                                }
                            }
                        }
                        else
                        {
                            /* ------------------------------------------------------ */
                            /* scatter kth column of the block into workspace X */
                            /* ------------------------------------------------------ */
                            oldcol = Q[k + k1];
                            pend = Ap[oldcol + 1];
                            for (p = Ap[oldcol]; p < pend; p++)
                            {
                                oldrow = Ai[p];
                                newrow = Pinv[oldrow] - k1;
                                if (newrow < 0 && poff < nzoff)
                                {
                                    /* entry in off-diagonal part */
                                    /* Offx [poff] = Az [p] / Rs [oldrow] */
                                    SCALE_DIV_ASSIGN(Offx[poff], Az[p], Rs[oldrow]);
                                    poff++;
                                }
                                else
                                {
                                    /* (newrow,k) is an entry in the block */
                                    /* X [newrow] = Az [p] / Rs [oldrow] */
                                    SCALE_DIV_ASSIGN(X[newrow], Az[p], Rs[oldrow]);
                                }
                            }

                            /* ------------------------------------------------------
                            */
                            /* compute kth column of U, and update kth column of
                            * A */
                            /* ------------------------------------------------------
                            */

                            GET_POINTER(LU, Uip, Ulen, Ui, Ux, k, ulen);
                            for (up = 0; up < ulen; up++)
                            {
                                j = Ui[up];
                                ujk = X[j];
                                /* X [j] = 0 */
                                CLEAR(X[j]);
                                Ux[up] = ujk;
                                GET_POINTER(LU, Lip, Llen, Li, Lx, j, llen);
                                for (p = 0; p < llen; p++)
                                {
                                    /* X [Li [p]] -= Lx [p] * ujk */
                                    MULT_SUB(X[Li[p]], Lx[p], ujk);
                                }
                            }
                            /* get the diagonal entry of U */
                            ukk = X[k];
                            /* X [k] = 0 */
                            CLEAR(X[k]);
                            /* singular case */
                            if (IS_ZERO(ukk))
                            {
                                /* matrix is numerically singular */
                                Common->status = KLU_SINGULAR;
                                if (Common->numerical_rank == EMPTY)
                                {
                                    Common->numerical_rank = k + k1;
                                    Common->singular_col = Q[k + k1];
                                }
                                if (Common->halt_if_singular)
                                {
                                    /* do not continue the factorization */
                                    return (FALSE);
                                }
                            }
                            /* pivot vadility testing */
                            else if (SCALAR_ABS(ukk) < Common->pivot_tol_fail)
                            {
                                /* pivot is too small */
                                Common->status = KLU_PIVOT_FAULT;
                                if (Common->halt_if_pivot_fails)
                                {
                                    /* do not continue the factorization */
                                    return (FALSE);
                                }
                            }
                            Udiag[k + k1] = ukk;
                            /* gather and divide by pivot to get kth column of L
                            */
                            GET_POINTER(LU, Lip, Llen, Li, Lx, k, llen);
                            for (p = 0; p < llen; p++)
                            {
                                i = Li[p];
                                DIV(Lx[p], X[i], ukk);
                                CLEAR(X[i]);
                            }
                        }
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* permute scale factors Rs according to pivotal row order */
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {
        for (k = 0; k < n; k++)
        {
            REAL(X[k]) = Rs[Pnum[k]];
        }
        for (k = 0; k < n; k++)
        {
            Rs[k] = REAL(X[k]);
        }
    }

#ifndef NDEBUG
    ASSERT(Numeric->Offp[n] == poff);
    ASSERT(Symbolic->nzoff == poff);
    PRINTF(("\n------------------- Off diagonal entries, new:\n"));
    ASSERT(KLU_valid(n, Numeric->Offp, Numeric->Offi, Offx));
    if (Common->status == KLU_OK)
    {
        PRINTF(("\n ########### KLU_BTF_REFACTOR done, nblocks %d\n", nblocks));
        for (block = 0; block < nblocks; block++)
        {
            k1 = R[block];
            k2 = R[block + 1];
            nk = k2 - k1;
            PRINTF(("\n================KLU_refactor output: k1 %d k2 %d nk %d\n", k1, k2, nk));
            if (nk == 1)
            {
                PRINTF(("singleton  "));
                PRINT_ENTRY(Udiag[k1]);
            }
            else
            {
                Lip = Numeric->Lip + k1;
                Llen = Numeric->Llen + k1;
                LU = (Unit *)Numeric->LUbx[block];
                PRINTF(("\n---- L block %d\n", block));
                ASSERT(KLU_valid_LU(nk, TRUE, Lip, Llen, LU));
                Uip = Numeric->Uip + k1;
                Ulen = Numeric->Ulen + k1;
                PRINTF(("\n---- U block %d\n", block));
                ASSERT(KLU_valid_LU(nk, FALSE, Uip, Ulen, LU));
            }
        }
    }
#endif
#ifdef KLU_PRINT
    static int counter = 0;
    if(counter == 0 || counter == 1000)
    {
        int n = Symbolic->n;
        int lnz = Numeric->lnz;
        int unz = Numeric->unz;
        int nzoff = Numeric->nzoff;
        int nb = Symbolic->nblocks;
        int *Lp, *Li, *Up, *Ui, *Fi, *Fp;
        double *Lx, *Ux, *Fx;
        int *P, *Q, *R;
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
        klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common);
        dumpKAll(Lx, Li, Lp, Ux, Ui, Up, Fx, Fi, Fp, P, Q, Numeric->path, Numeric->block_path, lnz, unz, n, nzoff, nb);
        dumpKA(Ax, Ai, Ap, n);
    }
    counter++;
#endif
    return (TRUE);
}
