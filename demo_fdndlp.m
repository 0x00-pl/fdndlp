
clc;
clear;
close all;

% *****************************************************
% Set path
% *****************************************************

addpath(genpath('_lib_'));
output_dir = 'wav_out/';
if ~exist(output_dir, 'dir')
   mkdir(output_dir); 
   disp(['mkdir ', output_dir])
end

%******************************************************
% Input and Output Configurations 
%******************************************************

filepath = 'wav_sample/';
sample_name = 'r_female_07_4ch.wav';
file_name = [filepath, sample_name];
out_name = [output_dir, strrep(sample_name, 'r_', 'd1_')];

sig_multi_mode = 1;
% sig_multi_mode = 0;

cfg.num_mic = 3;                            % the number of channels
cfg.num_out = 3;


%******************************************************
% Set Parameters
%******************************************************
cfg.K = 512;                               % the number of subbands
cfg.F = 2;                                 % over-sampling rate
% cfg.F = 1;                                 % over-sampling rate
cfg.N = cfg.K / cfg.F;                     % decimation factor


cfg.D1 = 2;                                % subband preditction delay                     
cfg.Lc1 = 30;                              % subband prediction order 
cfg.eps = 1e-4;                            % lower bound of rho(Normalizaton factor)
cfg.iterations = 2;



%******************************************************
% Read Audio Files and Processing
%******************************************************

disp('Read Audio Files...')

if sig_multi_mode
    disp(file_name)
    % audioread����
    % x������Ϊ m��n �������� m �Ƕ�ȡ����Ƶ��������n ���ļ��е���Ƶͨ����
    % fs��������
    [x, fs] = audioread(file_name);
    
    fid = fopen('wav_file_in.txt','w');
    fprintf(fid,'%e\n',x);
    fclose(fid);
    
    x = x(:,1:cfg.num_mic);
else
    x = [];
    for m = 1 : cfg.num_mic
        disp(file_name)
        filename1 = strrep(file_name, 'ch1', ['ch',num2str(m)]);
        [s, fs] = audioread(filename1);
        x = [x, s];
    end
end
cfg.fs = fs;

disp('Procssing...')
y = fdndlp(x, cfg);


% ***************************************************
% Output
% ***************************************************
disp(['write to file:',32, out_name])
udf.fig(x(:,1), fs);
udf.fig(y(:,1), fs);
% udf.play(x(:,1), fs);
% udf.play(y(:,1), fs);
audiowrite(out_name, y/max(max(abs(y))), fs);

rmpath(genpath('_lib_'));

