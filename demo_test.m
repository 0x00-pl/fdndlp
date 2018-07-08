
clc;
clear;
close all;

% Set path
addpath(genpath('_lib_'));
output_dir = 'wav_out/';
if ~exist(output_dir, 'dir')
   mkdir(output_dir); 
   disp(['mkdir ', output_dir])
end

% Input and Output Configurations 
filepath = 'wav_sample/';
sample_name = 'r_female_07_4ch.wav';
file_name = [filepath, sample_name];
out_name = [output_dir, strrep(sample_name, 'r_', 'd1_')];

% configs
cfg.fft_size = 512;
cfg.shift = cfg.fft_size/2;
cfg.fwindow = @hann;

% read audio
[ain, fs] = audioread(file_name);
cfg.audio_size_origin = size(ain,1);
cfg.N = floor((cfg.audio_size_origin-cfg.fft_size)/cfg.shift)+1;

fin = zeros(cfg.fft_size, cfg.N, size(ain, 2));
for i = 0:cfg.N-1
    tmp = ain(i*cfg.shift+1:(i*cfg.shift)+cfg.fft_size, :);
    ftmp = fft(tmp, cfg.fft_size, 1);
    fin(:, i+1, :) = reshape(ftmp, cfg.fft_size, 1, []);
end

fin1 = reshape(fin(:,:,1), cfg.fft_size, cfg.N);

imagesc(log(abs(real(reshape(fin(:,:,1), cfg.fft_size, cfg.N)))))
figure;
imagesc(log(abs(imag(reshape(fin(:,:,1), cfg.fft_size, cfg.N)))))
figure;
imagesc(log(abs(reshape(fin(:,:,1), cfg.fft_size, cfg.N))))
figure;
imagesc(real(log(reshape(fin(:,:,1), cfg.fft_size, cfg.N))))
figure;
imagesc(imag(log(reshape(fin(:,:,1), cfg.fft_size, cfg.N))))
figure;
imagesc(abs(log(reshape(fin(:,:,1), cfg.fft_size, cfg.N))))
figure;
plot(angle(mean(fin1)))

% fin = fft(ain, 512, 1);
% plot(abs(fin))
