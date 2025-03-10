//------------------------------------------------------------------------------
// SPEX/MATLAB/spex_cholesky_mex_soln: Use SPEX Chol within MATLAB
//------------------------------------------------------------------------------

// SPEX: (c) 2022-2024, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: The .c file defining the SPEX Chol MATLAB interfacee
 * This function defines: x = spex_cholesky_mex_soln (A, b, option)
 */


#include "SPEX_mex.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    //--------------------------------------------------------------------------
    // Initialize SPEX Chol library environment
    //--------------------------------------------------------------------------

    SPEX_info status ;
    SPEX_MEX_OK(SPEX_initialize_expert(mxMalloc, mxCalloc, mxRealloc, mxFree));
    SuiteSparse_config_printf_func_set (mexPrintf);

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    if (nargout > 1 || nargin < 2 || nargin > 3)
    {
        spex_mex_error (1, "Usage: x = SPEX_chol_soln (A,b,option)");
    }

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (mxIsComplex (pargin [0]) || mxIsComplex (pargin [1]))
    {
        spex_mex_error (1, "inputs must be real");
    }
    /***/
    if (!mxIsSparse (pargin [0]))     // Is the matrix sparse?
    {
        spex_mex_error (1, "first input must be sparse");
    }
    /**/
    if (mxIsSparse (pargin [1]))         // Is b sparse?
    {
        spex_mex_error (1, "second input must be full");
    }

    //--------------------------------------------------------------------------
    // get the input options
    //--------------------------------------------------------------------------

    SPEX_options option = NULL;
    SPEX_create_default_options(&option);
    if (option == NULL)
    {
        spex_mex_error (SPEX_OUT_OF_MEMORY, "");
    }

    spex_mex_options mexoptions ;
    if (nargin > 2) spex_mex_get_matlab_options (option, &mexoptions, pargin [2]);

    //--------------------------------------------------------------------------
    // get A and b
    //--------------------------------------------------------------------------

    SPEX_matrix A = NULL ;
    SPEX_matrix b = NULL ;
    spex_mex_get_A_and_b (&A, &b, pargin, option);
/**/
    if (option->print_level > 0)
    {
        printf ("\nScaled integer input matrix A:\n");
        SPEX_matrix_check (A, option);
    }

    if (option->print_level > 0)
    {
        printf ("\nScaled integer right-hand-side b:\n");
        SPEX_matrix_check (b, option);
    }

    //--------------------------------------------------------------------------
    // x = A\b via SPEX_Chol, returning result as SPEX_MPQ
    //--------------------------------------------------------------------------

    SPEX_matrix x = NULL ;
    SPEX_MEX_OK (SPEX_cholesky_backslash (&x, SPEX_MPQ, A, b, option));

    //--------------------------------------------------------------------------
    // print the result, if requested
    //--------------------------------------------------------------------------

    if (option->print_level > 0)
    {
        printf ("\nSolution x:\n");
        SPEX_matrix_check (x, option);
    }

    //--------------------------------------------------------------------------
    // return x to MATLAB
    //--------------------------------------------------------------------------

    if (mexoptions.solution == SPEX_SOLUTION_DOUBLE)
    {

        //----------------------------------------------------------------------
        // return x as a double MATLAB matrix
        //----------------------------------------------------------------------

        // convert x to double
        SPEX_matrix t = NULL ;
        SPEX_MEX_OK (SPEX_matrix_copy (&t, SPEX_DENSE, SPEX_FP64, x, option));
        SPEX_matrix_free (&x, NULL);
        x = t ;
        t = NULL ;

        // create an empty 0-by-0 MATLAB matrix and free its contents
        pargout [0] = mxCreateDoubleMatrix (0, 0, mxREAL);
        mxFree (mxGetDoubles (pargout [0]));

        // transplant x into the new MATLAB matrix and set its size
        mxSetDoubles (pargout [0], x->x.fp64);
        mxSetM (pargout [0], x->m);
        mxSetN (pargout [0], x->n);
        x->x.fp64 = NULL ;  // set to NULL so it is not freed by SPEX_matrix_free

    }
    else
    {

        //----------------------------------------------------------------------
        // return x as a cell array of strings
        //----------------------------------------------------------------------

        pargout [0] = mxCreateCellMatrix (x->m, x->n);
        int64_t mn = (x->m) * (x->n);
        for (int64_t p = 0 ; p < mn ; p++)
        {
            // convert x (i,j) into a C string
            char *s ;
            status = SPEX_mpfr_asprintf (&s, "%Qd", x->x.mpq [p]);
            if (status < 0)
            {
                spex_mex_error (1, "error converting x to string");
            }
            // convert the string into a MATLAB string and store in x {i,j}
            mxSetCell (pargout [0], p, mxCreateString (s));
            // free the C string
            SPEX_mpfr_free_str (s);
        }
    }

    //--------------------------------------------------------------------------
    // free workspace
    //--------------------------------------------------------------------------

    SPEX_matrix_free (&x, option);
    SPEX_matrix_free (&b, option);
    SPEX_matrix_free (&A, option);
    SPEX_FREE (option);
    SPEX_finalize ( );
}

#if defined ( __clang__ ) && defined ( CLANG_NEEDS_MAIN )
// when using clang inside MATLAB, it seems to require a "_main" symbol.
int main (void) { return (0) ; }
#endif

