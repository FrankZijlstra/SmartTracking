% ir_mex_build_table
% run matlab's "mex" command to "compile" the table interpolation code
% into mex files.
% only users on unsupported systems, e.g., PCs, will need to do this

% dir_current = pwd;
% dir_nufft = path_find_dir('nufft');
% dir_table = [dir_nufft filesep 'table'];
% cd(dir_table)

files = dir('*.c');
for n=1:length(files);
    movefile(files(n).name, [files(n).name(1:end-2) '.cpp']);
end

% mex compile command with optimization and c99 flag:
fun = @(f1, f2) mex('-O', 'CFLAGS="\$CFLAGS -std=c99"', f1, f2)

fun('interp1_table_adj_mex.cpp',	'interp1_table1_adj.cpp')
fun('interp1_table_mex.cpp',	'interp1_table1_for.cpp')
fun('interp2_table_adj_mex.cpp',	'interp2_table1_adj.cpp')
fun('interp2_table_mex.cpp',	'interp2_table1_for.cpp')
fun('interp3_table_adj_mex.cpp',	'interp3_table1_adj.cpp')
fun('interp3_table_mex.cpp',	'interp3_table1_for.cpp')

%{ old way
mex interp1_table_adj_mex.cpp	interp1_table1_adj.cpp
mex interp1_table_mex.cpp	interp1_table1_for.cpp
mex interp2_table_adj_mex.cpp	interp2_table1_adj.cpp
mex interp2_table_mex.cpp	interp2_table1_for.cpp
mex interp3_table_adj_mex.cpp	interp3_table1_adj.cpp
mex interp3_table_mex.cpp	interp3_table1_for.cpp
%}

cd(dir_current)
