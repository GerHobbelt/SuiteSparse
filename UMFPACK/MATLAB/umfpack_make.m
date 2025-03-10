function umfpack_make
%UMFPACK_MAKE to compile umfpack for use in MATLAB
%
% Compiles the umfpack mexFunction and then runs a simple demo.
%
% Example:
%   umfpack_make
%
% UMFPACK relies on AMD and its own built-in version of COLAMD for its ordering
% options.  The default is for UMFPACK to also use CHOLMOD, CCOLAMD, CAMD, and
% METIS for more ordering options as well.  This results in lower fill-in and
% higher performance.  METIS 5.1.0 is in ../../CHOLMOD/SuiteSparse_metis.
% METIS is optional; if not present, it is not used.
%
% See also: umfpack, umfpack_details, umfpack_report, umfpack_demo,
% and umfpack_simple.

% UMFPACK, Copyright (c) 2005-2022, Timothy A. Davis, All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0+

metis_path = '../../CHOLMOD/SuiteSparse_metis' ;
with_cholmod = exist (metis_path, 'dir') ;

details = 0 ;   % set to 1 to print out each mex command as it's executed

flags = '' ;
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    flags = ' -largeArrayDims' ;
end

% MATLAB 8.3.0 now has a -silent option to keep 'mex' from burbling too much
if (~verLessThan ('matlab', '8.3.0'))
    flags = ['-silent ' flags] ;
end

if (ispc)
    % MSVC does not define ssize_t
    flags = [flags ' -DNO_SSIZE_T'] ;
    % disable the SuiteSparse_config timer
    flags = ['-DNTIMER ' flags] ;
end

v = version ;

fprintf ('Compiling UMFPACK for MATLAB Version %s\n', v) ;

if (ispc)
    obj = 'obj' ;
else
    obj = 'o' ;
end

kk = 0 ;

%-------------------------------------------------------------------------------
% BLAS option
%-------------------------------------------------------------------------------

% This is exceedingly ugly.  The MATLAB mex command needs to be told where to
% find the LAPACK and BLAS libraries, which is a real portability nightmare.

if (ispc)
    % BLAS/LAPACK functions have no underscore on Windows
    flags = [flags ' -DBLAS_NO_UNDERSCORE'] ;
    if (verLessThan ('matlab', '7.5'))
        lapack = 'libmwlapack.lib' ;
    elseif (verLessThan ('matlab', '9.5'))
        lapack = 'libmwlapack.lib libmwblas.lib' ;
    else
        lapack = '-lmwlapack -lmwblas' ;
    end
else
    % BLAS/LAPACK functions have an underscore suffix
    flags = [flags ' -DBLAS_UNDERSCORE'] ;
    if (verLessThan ('matlab', '7.5'))
        lapack = '-lmwlapack' ;
    else
        lapack = '-lmwlapack -lmwblas' ;
    end
end

if (is64 && ~verLessThan ('matlab', '7.8'))
    % versions 7.8 and later on 64-bit platforms use a 64-bit BLAS
    fprintf ('with 64-bit BLAS\n') ;
    flags = [flags ' -DBLAS64'] ;
else
    % other versions of MATLAB use a 32-bit BLAS
    flags = [flags ' -DBLAS32'] ;
end

if (~(ispc || ismac))
    % for POSIX timing routine
    lapack = [lapack ' -lrt'] ;
end

%-------------------------------------------------------------------------------
% Source and include directories
%-------------------------------------------------------------------------------

umfdir = '../Source/' ;
amddir = '../../AMD/Source/' ;
incdir = ' -I. -I../Include -I../Source -I../../AMD/Include -I../../SuiteSparse_config' ;

if (with_cholmod)
    incdir = [incdir ' -I../../CCOLAMD/Include -I../../CAMD/Include ' ...
    ' -I../../CHOLMOD/Include -I../../COLAMD/Include -I../../CHOLMOD'] ;
    incdir = [incdir ' -I' metis_path '/include'] ;
    incdir = [incdir ' -I' metis_path '/GKlib'] ;
    incdir = [incdir ' -I' metis_path '/libmetis'] ;
end

%-------------------------------------------------------------------------------
% METIS options
%-------------------------------------------------------------------------------

if (with_cholmod)
    fprintf ('with CHOLMOD, CAMD, CCOLAMD, and METIS\n') ;
    flags = [' -DNSUPERNODAL -DNMODIFY -DNMATRIXOPS ' flags] ;
else
    fprintf ('without CHOLMOD, CAMD, CCOLAMD, and METIS\n') ;
    flags = [' -DNCHOLMOD ' flags] ;
end

%-------------------------------------------------------------------------------
% source files
%-------------------------------------------------------------------------------

% non-user-callable umf_*.[ch] files:
umfch = { 'assemble', 'blas3_update', ...
        'build_tuples', 'create_element', ...
        'dump', 'extend_front', 'garbage_collection', ...
        'get_memory', 'init_front', 'kernel', ...
        'kernel_init', 'kernel_wrapup', ...
        'local_search', 'lsolve', 'ltsolve', ...
        'mem_alloc_element', 'mem_alloc_head_block', ...
        'mem_alloc_tail_block', 'mem_free_tail_block', ...
        'mem_init_memoryspace', ...
        'report_vector', 'row_search', 'scale_column', ...
        'set_stats', 'solve', 'symbolic_usage', 'transpose', ...
        'tuple_lengths', 'usolve', 'utsolve', 'valid_numeric', ...
        'valid_symbolic', 'grow_front', 'start_front', ...
	'store_lu', 'scale' } ;

% non-user-callable umf_*.[ch] files, int versions only (no real/complex):
umfint = { 'analyze', 'apply_order', 'colamd', 'free', 'fsize', ...
        'is_permutation', 'malloc', 'realloc', 'report_perm', ...
	'singletons', 'cholmod' } ;

% non-user-callable and user-callable amd_*.[ch] files (int versions only):
amdsrc = { 'aat', '1', '2', 'dump', 'postorder', 'post_tree', 'defaults', ...
        'order', 'control', 'info', 'valid', 'preprocess' } ;

% user-callable umfpack_*.[ch] files (real/complex):
user = { 'col_to_triplet', 'defaults', 'free_numeric', ...
        'free_symbolic', 'get_numeric', 'get_lunz', ...
        'get_symbolic', 'get_determinant', 'numeric', 'qsymbolic', ...
        'report_control', 'report_info', 'report_matrix', ...
        'report_numeric', 'report_perm', 'report_status', ...
        'report_symbolic', 'report_triplet', ...
        'report_vector', 'solve', 'symbolic', ...
        'transpose', 'triplet_to_col', 'scale' ...
	'load_numeric', 'save_numeric', 'load_symbolic', 'save_symbolic' } ;

% user-callable umfpack_*.[ch], only one version
generic = { 'timer', 'tictoc' } ;

M = cell (0) ;

% add the SuiteSparse_time function
other_source = { '../../SuiteSparse_config/SuiteSparse_config' } ;

% add CHOLMOD and its supporting libraries
if (with_cholmod)

    ordering_src = { ...
        '../../CAMD/Source/camd_1', ...
        '../../CAMD/Source/camd_2', ...
        '../../CAMD/Source/camd_aat', ...
        '../../CAMD/Source/camd_control', ...
        '../../CAMD/Source/camd_defaults', ...
        '../../CAMD/Source/camd_dump', ...
        '../../CAMD/Source/camd_info', ...
        '../../CAMD/Source/camd_order', ...
        '../../CAMD/Source/camd_postorder', ...
        '../../CAMD/Source/camd_preprocess', ...
        '../../CAMD/Source/camd_valid', ...
        '../../COLAMD/Source/colamd', ...
        '../../CCOLAMD/Source/ccolamd' } ;

    cholmod_src = {
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
        '../../CHOLMOD/Cholesky/cholmod_l_amd', ...
        '../../CHOLMOD/Cholesky/cholmod_l_analyze', ...
        '../../CHOLMOD/Cholesky/cholmod_l_colamd', ...
        '../../CHOLMOD/Cholesky/cholmod_l_etree', ...
        '../../CHOLMOD/Cholesky/cholmod_l_postorder', ...
        '../../CHOLMOD/Cholesky/cholmod_l_rowcolcounts', ...
        '../../CHOLMOD/Partition/cholmod_l_ccolamd', ...
        '../../CHOLMOD/Partition/cholmod_l_csymamd', ...
        '../../CHOLMOD/Partition/cholmod_l_camd', ...
        '../../CHOLMOD/Partition/cholmod_l_metis', ...
        '../../CHOLMOD/Partition/cholmod_metis_wrapper', ...
        '../../CHOLMOD/Partition/cholmod_l_nesdis' } ;

    other_source = [other_source cholmod_src ordering_src] ;
end

%-------------------------------------------------------------------------------
% mex command
%-------------------------------------------------------------------------------

% with optimization:
mx = sprintf ('mex -O%s%s ', incdir, flags) ;
% no optimization:
% mx = sprintf ('mex -g %s%s%s ', incdir, flags) ;

%-------------------------------------------------------------------------------
% CHOLMOD, CAMD, C*OLAMD, METIS, SuiteSparse_config, and rand48 for Windows
%-------------------------------------------------------------------------------

for k = 1:length(other_source)
    ff = other_source {k} ;
    slash = strfind (ff, '/') ;
    slash = slash (end) + 1 ;
    o = ff (slash:end) ;
    kk = cmd (sprintf ('%s -DDLONG -c %s.c', mx, ff), kk, details) ;
    M {end+1} = [o '.' obj] ;
end

%-------------------------------------------------------------------------------
% Create the umfpack and amd2 mexFunctions for MATLAB (int versions only)
%-------------------------------------------------------------------------------

for k = 1:length(umfint)
    [M, kk] = make (M, '%s -DDLONG -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s.%s', mx, umfint {k}, umfint {k}, 'm', obj, umfdir, ...
	kk, details) ;
end

rules = { [mx ' -DDLONG'] , [mx ' -DZLONG'] } ;
kinds = { 'md', 'mz' } ;

for what = 1:2

    rule = rules {what} ;
    kind = kinds {what} ;

    [M, kk] = make (M, '%s -DCONJUGATE_SOLVE -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s.%s', rule, 'ltsolve', 'lhsolve', kind, obj, umfdir, ...
	kk, details) ;

    [M, kk] = make (M, '%s -DCONJUGATE_SOLVE -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s.%s', rule, 'utsolve', 'uhsolve', kind, obj, umfdir, ...
	kk, details) ;

    [M, kk] = make (M, '%s -DDO_MAP -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_map_nox.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDO_VALUES -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_nomap_x.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -c %sumf_%s.c', 'umf_%s.%s',  ...
        'umf_%s_%s_nomap_nox.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDO_MAP -DDO_VALUES -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_map_x.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DFIXQ -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s_fixq.%s', rule, 'assemble', 'assemble', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDROP -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s_drop.%s', rule, 'store_lu', 'store_lu', kind, obj, ...
	umfdir, kk, details) ;

    for k = 1:length(umfch)
        [M, kk] = make (M, '%s -c %sumf_%s.c', 'umf_%s.%s', 'umf_%s_%s.%s', ...
            rule, umfch {k}, umfch {k}, kind, obj, umfdir, kk, details) ;
    end

    [M, kk] = make (M, '%s -DWSOLVE -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
        'umfpack_%s_w%s.%s', rule, 'solve', 'solve', kind, obj, umfdir, ...
	kk, details) ;

    for k = 1:length(user)
        [M, kk] = make (M, '%s -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
            'umfpack_%s_%s.%s', rule, user {k}, user {k}, kind, obj, ...
	    umfdir, kk, details) ;
    end
end

for k = 1:length(generic)
    [M, kk] = make (M, '%s -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
	'umfpack_%s_%s.%s', mx, generic {k}, generic {k}, 'm', obj, ...
	umfdir, kk, details) ;
end

%----------------------------------------
% AMD routines (long only)
%----------------------------------------

for k = 1:length(amdsrc)
    [M, kk] = make (M, '%s -DDLONG -c %samd_%s.c', 'amd_%s.%s', ...
	'amd_%s_%s.%s', mx, amdsrc {k}, amdsrc {k}, 'm', obj, amddir, ...
	kk, details) ;
end

%----------------------------------------
% compile the umfpack mexFunction
%----------------------------------------

C = sprintf ('%s -output umfpack umfpackmex.c', mx) ;
for i = 1:length (M)
    C = [C ' ' (M {i})] ;   %#ok
end
C = [C ' ' lapack] ;
kk = cmd (C, kk, details) ;

%----------------------------------------
% delete the object files
%----------------------------------------

for i = 1:length (M)
    rmfile (M {i}) ;
end

%----------------------------------------
% compile the luflop mexFunction
%----------------------------------------

cmd (sprintf ('%s -output luflop luflopmex.c', mx), kk, details) ;

fprintf ('\nUMFPACK successfully compiled\n') ;

%===============================================================================
% end of umfpack_make
%===============================================================================


%-------------------------------------------------------------------------------

function rmfile (file)
% rmfile:  delete a file, but only if it exists
if (length (dir (file)) > 0)						    %#ok
    delete (file) ;
end

%-------------------------------------------------------------------------------

function cpfile (src, dst)
% cpfile:  copy the src file to the filename dst, overwriting dst if it exists
rmfile (dst)
if (length (dir (src)) == 0)	%#ok
    fprintf ('File does not exist: %s\n', src) ;
    error ('File does not exist') ;
end
try
    copyfile (src, dst) ;
catch ME
    % ignore errors of the form "cp: preserving permissions: ...
    % Operation not supported".  rethrow all other errors.
    if (isempty (strfind (ME.message, 'Operation not supported')))
        rethrow (ME) ;
    end
end

%-------------------------------------------------------------------------------

function mvfile (src, dst)
% mvfile:  move the src file to the filename dst, overwriting dst if it exists
cpfile (src, dst) ;
rmfile (src) ;

%-------------------------------------------------------------------------------

function kk = cmd (s, kk, details)
%CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
	fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;

%-------------------------------------------------------------------------------

function [M, kk] = make (M, s, src, dst, rule, file1, file2, kind, obj, ...
    srcdir, kk, details)
% make:  execute a "make" command for a source file
kk = cmd (sprintf (s, rule, srcdir, file1), kk, details) ;
src = sprintf (src, file1, obj) ;
dst = sprintf (dst, kind, file2, obj) ;
mvfile (src, dst) ;
M {end + 1} = dst ;
