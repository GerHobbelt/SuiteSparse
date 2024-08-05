function paru_make (try_intel)
%PARU_MAKE compiles the ParU mexFunction
%
% Usage: paru_make
%
% Example:
%
%   paru_make
%   load west0479
%   A = west0479 ;
%   n = size (A,1) ;
%   b = rand (n,1) ;
%   x1 = A\b ;
%   norm (A*x1-b)
%   x2 = paru (A,b) ;
%   norm (A*x2-b)
%
% You must type the paru_make command while in the ParU/MATLAB directory.
%
% For best performance, paru relies on functions unique to the Intel MKL BLAS.
% An optional input, paru_make(try_intel), is true by default.  paru_make
% detects the BLAS library used by MATLAB and then attempts to use functions
% unique to the Intel MKL BLAS library (mkl*set_num_threads_local).  This may
% fail when paru is compiled (in which case paru is compiled with try_intel
% false).  If paru fails when it runs, with a link error reporting that an an
% mkl_* routine is not found, use paru_make(false) to disable the Intel MKL
% functions.
%
% See also paru, paru_demo, paru_many, paru_tiny, mldivide.

% ParU, Copyright (c) 2022-2024, Mohsen Aznaveh and Timothy A. Davis,
% All Rights Reserved.
% SPDX-License-Identifier: GPL-3.0-or-later

if (nargin < 1)
    try_intel = (~ismac && isunix) ;
end
if (try_intel)
    v = version ('-blas') ;
    try_intel = contains (v, 'Intel') ;
end
if (try_intel)
    try
        paru_make_helper (true) ;
    catch
        fprintf ('\nCompilation with Intel MKL specific methods failed.\n') ;
        fprintf ('Retrying with standard BLAS API.\n') ;
        paru_make_helper (false) ;
    end
else
    paru_make_helper (false) ;
end
end

%-------------------------------------------------------------------------------

function paru_make_helper (try_intel)
fprintf ('Compiling ParU for MATLAB.\n') ;
if (try_intel)
    fprintf ('Using the Intel MKL BLAS\n') ;
end

% -R2018a: interleaved complex is required
flags = '-O -R2018a -silent ' ;

% developer-only flag, for testing only, not intended for end users:
flags = [flags ' -DDEVELOPER=1 ' ] ;

if (ispc)
    % MSVC does not define ssize_t
    flags = [flags ' -DNO_SSIZE_T -DNTIMER -DBLAS_NO_UNDERSCORE '] ;
end

if (ismac)
    flags = [flags ' -DCLANG_NEEDS_MAIN'] ;
end

libs = '-lmwlapack -lmwblas' ;
if (~(ispc || ismac))
    % for POSIX timing routine
    libs = [libs ' -lrt'] ;
end

openmp = '' ;
if (~ismac && isunix)
    openmp = ' -fopenmp' ;
end
flags = sprintf ('%s CFLAGS=''-fPIC %s'' LDFLAGS=''%s''', flags, openmp, openmp) ;

if (ispc)
    obj = 'obj' ;
else
    obj = 'o' ;
end

% check if MATLAB has the MKL Intel BLAS on Linux or Mac
if (try_intel && (ismac || isunix))
    v = version ('-blas') ;
    if (contains (v, 'Intel'))
        flags = [flags ' -DBLAS_Intel10_64ilp'] ;
    end
end

% use all of CHOLMOD except for the Modify module
flags = [flags ' -DNMODIFY -DBLAS64' ] ;
if (ispc)
    % MSVC does not define ssize_t
    flags = [flags ' -DNO_SSIZE_T'] ;
end

%-------------------------------------------------------------------------------

include = '-I../../include -I../Include' ;
include = [include ' -I../../SuiteSparse_config '] ;
include = [include ' -I../../AMD/Source  -I../../AMD/Include  '] ;
include = [include ' -I../../CAMD/Source -I../../CAMD/Include '] ;
include = [include ' -I../../COLAMD/Source  -I../../COLAMD/Include  '] ;
include = [include ' -I../../CCOLAMD/Source -I../../CCOLAMD/Include '] ;
include = [include ' -I../../UMFPACK/Source -I../../UMFPACK/Include '] ;
include = [include ' -I../../CHOLMOD/Include '] ;
include = [include ' -I../../CHOLMOD/Check '] ;
include = [include ' -I../../CHOLMOD/Cholesky '] ;
include = [include ' -I../../CHOLMOD/Utility '] ;
include = [include ' -I../../CHOLMOD/MATLAB '] ;
include = [include ' -I../../CHOLMOD/MatrixOps '] ;
include = [include ' -I../../CHOLMOD/Partition '] ;
include = [include ' -I../../CHOLMOD '] ;
include = [include ' -I../../CHOLMOD/SuiteSparse_metis/GKlib '] ;
include = [include ' -I../../CHOLMOD/SuiteSparse_metis/libmetis '] ;
include = [include ' -I../../CHOLMOD/SuiteSparse_metis/include '] ;
include = [include ' -I../../CHOLMOD/Supernodal '] ;

suitesparse_src = { ...
    '../../SuiteSparse_config/SuiteSparse_config', ...
    '../../AMD/Source/amd_l1', ...
    '../../AMD/Source/amd_l2', ...
    '../../AMD/Source/amd_l_aat', ...
    '../../AMD/Source/amd_l_control', ...
    '../../AMD/Source/amd_l_defaults', ...
    '../../AMD/Source/amd_l_dump', ...
    '../../AMD/Source/amd_l_info', ...
    '../../AMD/Source/amd_l_order', ...
    '../../AMD/Source/amd_l_postorder', ...
    '../../AMD/Source/amd_l_post_tree', ...
    '../../AMD/Source/amd_l_preprocess', ...
    '../../AMD/Source/amd_l_valid', ...
    '../../CAMD/Source/camd_l1', ...
    '../../CAMD/Source/camd_l2', ...
    '../../CAMD/Source/camd_l_aat', ...
    '../../CAMD/Source/camd_l_control', ...
    '../../CAMD/Source/camd_l_defaults', ...
    '../../CAMD/Source/camd_l_dump', ...
    '../../CAMD/Source/camd_l_info', ...
    '../../CAMD/Source/camd_l_order', ...
    '../../CAMD/Source/camd_l_postorder', ...
    '../../CAMD/Source/camd_l_preprocess', ...
    '../../CAMD/Source/camd_l_valid', ...
    '../../COLAMD/Source/colamd_l', ...
    '../../CCOLAMD/Source/ccolamd_l', ...
    '../../CHOLMOD/MATLAB/sputil2', ...
    '../../CHOLMOD/Utility/cholmod_l_aat', ...
    '../../CHOLMOD/Utility/cholmod_l_add', ...
    '../../CHOLMOD/Utility/cholmod_l_add_size_t', ...
    '../../CHOLMOD/Utility/cholmod_l_allocate_dense', ...
    '../../CHOLMOD/Utility/cholmod_l_allocate_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_allocate_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_allocate_triplet', ...
    '../../CHOLMOD/Utility/cholmod_l_allocate_work', ...
    '../../CHOLMOD/Utility/cholmod_l_alloc_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_alloc_work', ...
    '../../CHOLMOD/Utility/cholmod_l_band', ...
    '../../CHOLMOD/Utility/cholmod_l_band_nnz', ...
    '../../CHOLMOD/Utility/cholmod_l_calloc', ...
    '../../CHOLMOD/Utility/cholmod_l_change_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_clear_flag', ...
    '../../CHOLMOD/Utility/cholmod_l_copy', ...
    '../../CHOLMOD/Utility/cholmod_l_copy_dense2', ...
    '../../CHOLMOD/Utility/cholmod_l_copy_dense', ...
    '../../CHOLMOD/Utility/cholmod_l_copy_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_copy_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_copy_triplet', ...
    '../../CHOLMOD/Utility/cholmod_l_cumsum', ...
    '../../CHOLMOD/Utility/cholmod_l_dbound', ...
    '../../CHOLMOD/Utility/cholmod_l_defaults', ...
    '../../CHOLMOD/Utility/cholmod_l_dense_nnz', ...
    '../../CHOLMOD/Utility/cholmod_l_dense_to_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_divcomplex', ...
    '../../CHOLMOD/Utility/cholmod_l_ensure_dense', ...
    '../../CHOLMOD/Utility/cholmod_l_error', ...
    '../../CHOLMOD/Utility/cholmod_l_eye', ...
    '../../CHOLMOD/Utility/cholmod_l_factor_to_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_finish', ...
    '../../CHOLMOD/Utility/cholmod_l_free', ...
    '../../CHOLMOD/Utility/cholmod_l_free_dense', ...
    '../../CHOLMOD/Utility/cholmod_l_free_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_free_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_free_triplet', ...
    '../../CHOLMOD/Utility/cholmod_l_free_work', ...
    '../../CHOLMOD/Utility/cholmod_l_hypot', ...
    '../../CHOLMOD/Utility/cholmod_l_malloc', ...
    '../../CHOLMOD/Utility/cholmod_l_maxrank', ...
    '../../CHOLMOD/Utility/cholmod_l_mult_size_t', ...
    '../../CHOLMOD/Utility/cholmod_l_nnz', ...
    '../../CHOLMOD/Utility/cholmod_l_ones', ...
    '../../CHOLMOD/Utility/cholmod_l_pack_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_ptranspose', ...
    '../../CHOLMOD/Utility/cholmod_l_reallocate_column', ...
    '../../CHOLMOD/Utility/cholmod_l_reallocate_factor', ...
    '../../CHOLMOD/Utility/cholmod_l_reallocate_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_reallocate_triplet', ...
    '../../CHOLMOD/Utility/cholmod_l_realloc', ...
    '../../CHOLMOD/Utility/cholmod_l_realloc_multiple', ...
    '../../CHOLMOD/Utility/cholmod_l_sbound', ...
    '../../CHOLMOD/Utility/cholmod_l_score_comp', ...
    '../../CHOLMOD/Utility/cholmod_l_set_empty', ...
    '../../CHOLMOD/Utility/cholmod_l_sort', ...
    '../../CHOLMOD/Utility/cholmod_l_sparse_to_dense', ...
    '../../CHOLMOD/Utility/cholmod_l_sparse_to_triplet', ...
    '../../CHOLMOD/Utility/cholmod_l_speye', ...
    '../../CHOLMOD/Utility/cholmod_l_spzeros', ...
    '../../CHOLMOD/Utility/cholmod_l_start', ...
    '../../CHOLMOD/Utility/cholmod_l_transpose', ...
    '../../CHOLMOD/Utility/cholmod_l_transpose_sym', ...
    '../../CHOLMOD/Utility/cholmod_l_transpose_unsym', ...
    '../../CHOLMOD/Utility/cholmod_l_triplet_to_sparse', ...
    '../../CHOLMOD/Utility/cholmod_l_version', ...
    '../../CHOLMOD/Utility/cholmod_l_xtype', ...
    '../../CHOLMOD/Utility/cholmod_l_zeros', ...
    '../../CHOLMOD/Utility/cholmod_mult_uint64_t', ...
    '../../CHOLMOD/Utility/cholmod_memdebug', ...
    '../../CHOLMOD/Check/cholmod_l_check', ...
    '../../CHOLMOD/Check/cholmod_l_read', ...
    '../../CHOLMOD/Check/cholmod_l_write', ...
    '../../CHOLMOD/Cholesky/cholmod_l_amd', ...
    '../../CHOLMOD/Cholesky/cholmod_l_analyze', ...
    '../../CHOLMOD/Cholesky/cholmod_l_colamd', ...
    '../../CHOLMOD/Cholesky/cholmod_l_etree', ...
    '../../CHOLMOD/Cholesky/cholmod_l_factorize', ...
    '../../CHOLMOD/Cholesky/cholmod_l_postorder', ...
    '../../CHOLMOD/Cholesky/cholmod_l_rcond', ...
    '../../CHOLMOD/Cholesky/cholmod_l_resymbol', ...
    '../../CHOLMOD/Cholesky/cholmod_l_rowcolcounts', ...
    '../../CHOLMOD/Cholesky/cholmod_l_rowfac', ...
    '../../CHOLMOD/Cholesky/cholmod_l_solve', ...
    '../../CHOLMOD/Cholesky/cholmod_l_spsolve', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_drop', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_horzcat', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_norm', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_scale', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_sdmult', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_ssmult', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_submatrix', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_vertcat', ...
    '../../CHOLMOD/MatrixOps/cholmod_l_symmetry', ...
    '../../CHOLMOD/Supernodal/cholmod_l_super_numeric', ...
    '../../CHOLMOD/Supernodal/cholmod_l_super_solve', ...
    '../../CHOLMOD/Supernodal/cholmod_l_super_symbolic', ...
    '../../CHOLMOD/Partition/cholmod_metis_wrapper', ...
    '../../CHOLMOD/Partition/cholmod_l_ccolamd', ...
    '../../CHOLMOD/Partition/cholmod_l_csymamd', ...
    '../../CHOLMOD/Partition/cholmod_l_camd', ...
    '../../CHOLMOD/Partition/cholmod_l_metis', ...
    '../../CHOLMOD/Partition/cholmod_l_nesdis', ...
    '../../UMFPACK/Source2/umf_l_analyze', ...
    '../../UMFPACK/Source2/umf_l_apply_order', ...
    '../../UMFPACK/Source2/umf_l_cholmod', ...
    '../../UMFPACK/Source2/umf_l_colamd', ...
    '../../UMFPACK/Source2/umf_l_free', ...
    '../../UMFPACK/Source2/umf_l_fsize', ...
    '../../UMFPACK/Source2/umf_l_is_permutation', ...
    '../../UMFPACK/Source2/umf_l_malloc', ...
    '../../UMFPACK/Source2/umf_l_realloc', ...
    '../../UMFPACK/Source2/umf_l_report_perm', ...
    '../../UMFPACK/Source2/umf_l_singletons', ...
    '../../UMFPACK/Source2/umf_dl_assemble', ...
    '../../UMFPACK/Source2/umf_dl_assemble_fixq', ...
    '../../UMFPACK/Source2/umf_dl_blas3_update', ...
    '../../UMFPACK/Source2/umf_dl_build_tuples', ...
    '../../UMFPACK/Source2/umf_dl_create_element', ...
    '../../UMFPACK/Source2/umf_dl_dump', ...
    '../../UMFPACK/Source2/umf_dl_extend_front', ...
    '../../UMFPACK/Source2/umf_dl_garbage_collection', ...
    '../../UMFPACK/Source2/umf_dl_get_memory', ...
    '../../UMFPACK/Source2/umf_dl_grow_front', ...
    '../../UMFPACK/Source2/umf_dl_init_front', ...
    '../../UMFPACK/Source2/umf_dl_kernel', ...
    '../../UMFPACK/Source2/umf_dl_kernel_init', ...
    '../../UMFPACK/Source2/umf_dl_kernel_wrapup', ...
    '../../UMFPACK/Source2/umf_dl_lhsolve', ...
    '../../UMFPACK/Source2/umf_dl_local_search', ...
    '../../UMFPACK/Source2/umf_dl_lsolve', ...
    '../../UMFPACK/Source2/umf_dl_ltsolve', ...
    '../../UMFPACK/Source2/umf_dl_mem_alloc_element', ...
    '../../UMFPACK/Source2/umf_dl_mem_alloc_head_block', ...
    '../../UMFPACK/Source2/umf_dl_mem_alloc_tail_block', ...
    '../../UMFPACK/Source2/umf_dl_mem_free_tail_block', ...
    '../../UMFPACK/Source2/umf_dl_mem_init_memoryspace', ...
    '../../UMFPACK/Source2/umf_dl_report_vector', ...
    '../../UMFPACK/Source2/umf_dl_row_search', ...
    '../../UMFPACK/Source2/umf_dl_scale', ...
    '../../UMFPACK/Source2/umf_dl_scale_column', ...
    '../../UMFPACK/Source2/umf_dl_set_stats', ...
    '../../UMFPACK/Source2/umf_dl_solve', ...
    '../../UMFPACK/Source2/umf_dl_start_front', ...
    '../../UMFPACK/Source2/umf_dl_store_lu', ...
    '../../UMFPACK/Source2/umf_dl_store_lu_drop', ...
    '../../UMFPACK/Source2/umf_dl_symbolic_usage', ...
    '../../UMFPACK/Source2/umf_dl_transpose', ...
    '../../UMFPACK/Source2/umf_dl_triplet_map_nox', ...
    '../../UMFPACK/Source2/umf_dl_triplet_map_x', ...
    '../../UMFPACK/Source2/umf_dl_triplet_nomap_nox', ...
    '../../UMFPACK/Source2/umf_dl_triplet_nomap_x', ...
    '../../UMFPACK/Source2/umf_dl_tuple_lengths', ...
    '../../UMFPACK/Source2/umf_dl_uhsolve', ...
    '../../UMFPACK/Source2/umf_dl_usolve', ...
    '../../UMFPACK/Source2/umf_dl_utsolve', ...
    '../../UMFPACK/Source2/umf_dl_valid_numeric', ...
    '../../UMFPACK/Source2/umf_dl_valid_symbolic', ...
    '../../UMFPACK/Source2/umfpack_gn_tictoc', ...
    '../../UMFPACK/Source2/umfpack_gn_timer', ...
    '../../UMFPACK/Source2/umfpack_dl_col_to_triplet', ...
    '../../UMFPACK/Source2/umfpack_dl_defaults', ...
    '../../UMFPACK/Source2/umfpack_dl_free_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_free_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_get_determinant', ...
    '../../UMFPACK/Source2/umfpack_dl_get_lunz', ...
    '../../UMFPACK/Source2/umfpack_dl_get_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_get_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_load_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_load_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_qsymbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_report_control', ...
    '../../UMFPACK/Source2/umfpack_dl_report_info', ...
    '../../UMFPACK/Source2/umfpack_dl_report_matrix', ...
    '../../UMFPACK/Source2/umfpack_dl_report_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_report_perm', ...
    '../../UMFPACK/Source2/umfpack_dl_report_status', ...
    '../../UMFPACK/Source2/umfpack_dl_report_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_report_triplet', ...
    '../../UMFPACK/Source2/umfpack_dl_report_vector', ...
    '../../UMFPACK/Source2/umfpack_dl_save_numeric', ...
    '../../UMFPACK/Source2/umfpack_dl_save_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_scale', ...
    '../../UMFPACK/Source2/umfpack_dl_solve', ...
    '../../UMFPACK/Source2/umfpack_dl_symbolic', ...
    '../../UMFPACK/Source2/umfpack_dl_transpose', ...
    '../../UMFPACK/Source2/umfpack_dl_triplet_to_col', ...
    '../../UMFPACK/Source2/umfpack_dl_wsolve' } ;

paru_src = {
    '../Source/ParU_Factorize', ...
    '../Source/ParU_Analyze', ...
    '../Source/paru_assemble', ...
    '../Source/paru_assemble_row2U', ...
    '../Source/paru_bin_search', ...
    '../Source/ParU_C', ...
    '../Source/paru_create_element', ...
    '../Source/paru_cumsum', ...
    '../Source/paru_dgemm', ...
    '../Source/paru_diag_update', ...
    '../Source/paru_exec_tasks', ...
    '../Source/paru_finalize_perm', ...
    '../Source/ParU_FreeNumeric', ...
    '../Source/ParU_FreeSymbolic', ...
    '../Source/ParU_FreeControl', ...
    '../Source/ParU_InitControl', ...
    '../Source/paru_front', ...
    '../Source/paru_fs_factorize', ...
    '../Source/paru_full_summed', ...
    '../Source/paru_gaxpy', ...
    '../Source/ParU_Get', ...
    '../Source/paru_hash', ...
    '../Source/paru_heap', ...
    '../Source/paru_init_rel', ...
    '../Source/paru_init_rowFronts', ...
    '../Source/paru_intersection', ...
    '../Source/ParU_InvPerm', ...
    '../Source/ParU_LSolve', ...
    '../Source/paru_mem', ...
    '../Source/paru_free_work', ...
    '../Source/paru_memcpy', ...
    '../Source/paru_memset', ...
    '../Source/paru_norms', ...
    '../Source/paru_nthreads', ...
    '../Source/ParU_Perm', ...
    '../Source/paru_pivotal', ...
    '../Source/paru_prior_assemble', ...
    '../Source/ParU_Residual', ...
    '../Source/ParU_Set', ...
    '../Source/ParU_Solve', ...
    '../Source/paru_tasked_dgemm', ...
    '../Source/paru_tasked_dtrsm', ...
    '../Source/paru_dtrsm', ...
    '../Source/paru_tuples', ...
    '../Source/paru_umfpack_info', ...
    '../Source/paru_update_rel_ind', ...
    '../Source/paru_update_rowDeg', ...
    '../Source/ParU_USolve', ...
    '../Source/ParU_Version' } ;

obj_files = ' ' ;
objs = { } ;

% compile each ParU C++ file
for k = 1:length (paru_src)
    src = paru_src {k} ;
    s = sprintf ('mex -c %s %s %s.cpp', flags, include, src) ;
    % fprintf ('%s\n', src) ;
    fprintf ('.') ;
    eval (s) ;
    slash = strfind (src, '/') ;
    slash = slash (end) + 1 ;
    o = src (slash:end) ;
    obj_files = [ obj_files ' ' o '.' obj ] ; %#ok<AGROW>
    objs {end+1} = [o '.' obj] ;
    if (mod (k,60) == 0)
        fprintf ('\n') ;
    end
end
fprintf ('\n') ;

% compile each SuiteSparse C file (SuiteSparse_config, AMD, COLAMD,
% CAMD, CCOLAMD, UMFPACK, CHOLMOD):
for k = 1:length (suitesparse_src)
    src = suitesparse_src {k} ;
    s = sprintf ('mex -c %s %s %s.c', flags, include, src) ;
    % fprintf ('%s\n', src) ;
    fprintf ('.') ;
    eval (s) ;
    slash = strfind (src, '/') ;
    slash = slash (end) + 1 ;
    o = src (slash:end) ;
    obj_files = [ obj_files ' ' o '.' obj ] ; %#ok<AGROW>
    objs {end+1} = [o '.' obj] ;
    if (mod (k,60) == 0)
        fprintf ('\n') ;
    end
end

% compile the paru mexFunction
fprintf (':') ;
s = sprintf ('mex %s -O %s paru.c %s %s', flags, include, obj_files, libs) ;
% fprintf ('%s\n', s) ;
eval (s) ;
fprintf ('\n') ;

% delete the object files
for k = 1:length (objs)
    o = objs {k} ;
    if (length (dir (o)) > 0)
        delete (o) ;
    end
end

% try a quick demo
paru_tiny
paru_demo

end
