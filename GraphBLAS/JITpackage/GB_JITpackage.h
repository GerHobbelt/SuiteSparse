//------------------------------------------------------------------------------
// GB_JITpackage.h: definitions to package GraphBLAS source code for the JIT
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2024, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//------------------------------------------------------------------------------

#include "GraphBLAS.h"

typedef struct
{
    size_t uncompressed_size ;  // size of the uncompressed file
    size_t compressed_size ;    // size of the compressed blob
    uint8_t *blob ;             // the compressed blob
    char *filename ;            // filename, including the suffix
}
GB_JITpackage_index_struct ;

GB_GLOBAL int GB_JITpackage_nfiles ;
GB_GLOBAL GB_JITpackage_index_struct GB_JITpackage_index [ ] ;

