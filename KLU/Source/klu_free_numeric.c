/* ========================================================================== */
/* === KLU_free_numeric ===================================================== */
/* ========================================================================== */

/* Free the KLU Numeric object. */

#include "klu_internal.h"

Int KLU_free_numeric
(
    KLU_numeric **NumericHandle,
    KLU_common  *Common
)
{
    KLU_numeric *Numeric ;
    Unit **LUbx ;
    size_t *LUsize ;
    Int block, n, nzoff, nblocks ;

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (NumericHandle == NULL || *NumericHandle == NULL)
    {
        return (TRUE) ;
    }

    Numeric = *NumericHandle ;

    n = Numeric->n ;
    nzoff = Numeric->nzoff ;
    nblocks = Numeric->nblocks ;
    LUsize = Numeric->LUsize ;

    LUbx = (Unit **) Numeric->LUbx ;
    if (LUbx != NULL)
    {
        for (block = 0 ; block < nblocks ; block++)
        {
            KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
                sizeof (Unit), Common) ;
        }
    }

    KLU_free (Numeric->Pnum, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Offp, n+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offi, nzoff+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offx, nzoff+1, sizeof (Entry), Common) ;

    KLU_free (Numeric->Lip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Llen, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Uip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Ulen, n, sizeof (Int), Common) ;

    KLU_free (Numeric->LUsize, nblocks, sizeof (size_t), Common) ;

    KLU_free (Numeric->LUbx, nblocks, sizeof (Unit *), Common) ;

    KLU_free (Numeric->Udiag, n, sizeof (Entry), Common) ;

    KLU_free (Numeric->Rs,   n, sizeof (double), Common) ;
    KLU_free (Numeric->Pinv, n, sizeof (Int), Common) ;

    KLU_free (Numeric->Work, Numeric->worksize, 1, Common) ;

    /* Numeric->path or Numeric->bpath only allocated if klu_compute_path
     * is called, i.e. partial refactorization is used */
    if(Numeric->path)
    {
        KLU_free (Numeric->path, n, sizeof (int), Common) ;
    }
    if(Numeric->bpath)
    {
        KLU_free (Numeric->bpath, Numeric->nblocks, sizeof (int), Common);
    }
    KLU_free (Numeric, 1, sizeof (KLU_numeric), Common) ;

    *NumericHandle = NULL ;
    return (TRUE) ;
}
