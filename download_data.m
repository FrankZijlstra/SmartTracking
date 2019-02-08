% Download all datasets that can be reconstructed with run_tracking.m

clear all
close all

baseUrl = 'https://www.isi.uu.nl/People/Frank/smart_data/';

dataNames = {'spheres_linear', ...
             'spheres_rotate', ...
             'spheres_chaos', ...
             'needle_straight', ...
             'needle_45degrees'};

if (~exist('./data', 'dir'))
    mkdir('./data');
end

for I=1:length(dataNames)
    fprintf('Downloading: %s...\n', dataNames{I});
    websave(['./data/data_' dataNames{I} '.mat'], [baseUrl 'data_' dataNames{I} '.mat']);
    websave(['./data/template_' dataNames{I} '.mat'], [baseUrl 'template_' dataNames{I} '.mat']);
end
