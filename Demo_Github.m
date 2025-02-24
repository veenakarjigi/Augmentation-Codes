clc;clear;close all

tic

%% loading the trained model
addpath("D:\VK_HMC_Spectrogram_Augmentation_Work\VK_HMC_AUGMENT_FEB_2025");
folderPath = "D:\VK_HMC_Spectrogram_Augmentation_Work\VK_HMC_AUGMENT_FEB_2025\demo\unknown";
load D:\VK_HMC_Spectrogram_Augmentation_Work\models\augmented_10.mat;

duration = 1; % Recording duration in seconds
recorder = audiorecorder(fs, 16, 1); % Create audio recorder object
disp('Start speaking now...');
recordblocking(recorder, duration); 
disp('Recording complete.')

%% generation of spectrograms

% (a) STFT Spectrograms
fs = 48000; % Sampling frequency (Adjust if needed)
audioData=audioread("D:\VK_HMC_Spectrogram_Augmentation_Work\VK_HMC_AUGMENT_FEB_2025\two_words\Two Words_data\Aane\F01_002_1.wav");
audioData(:,2)=[];
N = 4096; % Number of FFT point
win_size = [240 480 720 960 1200 1440 1680 1920];  


% (1) STFT Spectrogram 
for i = 1:length(win_size)
overlap = win_size(i)/2; 
[Y,F,T] = spectrogram(audioData,win_size(i),overlap,N,fs);
Y((N/2)+1,:)=[];
Y_magnitude = abs(Y);
STFT_spec = 20*log10(Y_magnitude);  
STFT_spec1 = imresize(STFT_spec,[1000,1000]);
imageNameSTFT = sprintf("demo_stft(%d).png",i);
fullFilePathSTFT = fullfile(folderPath, imageNameSTFT);
imagesc(flipud(STFT_spec1));
saveas(gca,fullFilePathSTFT)

 % (2) STFT Mel Spectrogram 
 STFT_mfb_array = [];
 [temp_array,linear_centre_freq] = mfb_temp_array_fn(fs,N);
 [R,C] = size(Y_magnitude);
 for k = 1:C
     Xm = Y_magnitude(:,k);   
     res = temp_array*Xm;
     STFT_mfb_array(:,k) = res;
 end
            
 STFT_mel= 20*log10(STFT_mfb_array);
 STFT_mel1= imresize(STFT_mel, [1000,1000]);
 imagesc(flipud(STFT_mel1));
 imageNameSTFTmel = sprintf("demo_stftmel(%d).png",i);
 fullFilePathSTFTmel = fullfile(folderPath, imageNameSTFTmel);
 saveas(gca,fullFilePathSTFTmel)

% (3) STFT Mel PE Spectrogram 
    
[y1] = postaud(STFT_mfb_array,fs/2,'bark');
log_y1 = 20*log10(y1);
log_y2 = imresize(log_y1, [1000,1000]);
imagesc(flipud(log_y2));
imageNameSTFTmelpe = sprintf("demo_stftmelpe(%d).png",i);
fullFilePathSTFTmelpe = fullfile(folderPath, imageNameSTFTmelpe);
saveas(gca,fullFilePathSTFTmelpe)
end
        
% (b) SFF Spectrogram-Compute SFF Spectrum and Get the SFF spectrum only till fs/2 (nyquist frequency)
nfft = N;                               % number of FFT points
step_Hz =[20 80 200 800 1000 1200 1800 2000];        % Note if fhz/hfhz are used this is redundant
pflag = 0;                              % plot flag

for q = 1:length(step_Hz)
    nfftby2 = round(nfft/2);                % half-spectrum
    fhz = linspace(0,fs,nfft);              % full spectrum in Hz
    hfhz = linspace(0,fs/2,nfftby2);        % half-spectrum in Hz

    [sffspec] = singleFrequencyFilt(audioData,fs,step_Hz(q),0,hfhz);
    sffspec(sffspec <= 1e-9) = 1e-9;
    magnitude_sffspec = abs(sffspec);
    log_sffspec = 20*log10(magnitude_sffspec);   % Compute Log-Spectrum
    log_sffspec1 = imresize(log_sffspec, [1000,1000]);
    imagesc((log_sffspec1));
    imageNameSFF = sprintf("demo_sff(%d).png",q);
    fullFilePathsff = fullfile(folderPath, imageNameSFF);    
    saveas(gca,fullFilePathsff)
end

  %% (c) Constant Q Transform

        bins = [2 4 8 16 24 32 64 96];% CQT resolution [bins/octave]
        for q = 1:length(bins)
            fmax = fs/2;                    % highest frequency analyzed
            fmin = 100;                     % lowest frequency of interest (CQT bins will start immediatly above fmin)

            [cfs,f,g,fshifts,fintervals,bw] = cqt(audioData,'SamplingFrequency',fs,'BinsPerOctave',bins(q),'FrequencyLimits',[fmin fmax],'Window','hamming');
            cfs = cfs.c;
            C_Q_T = 20*log10(abs(cfs(1:length(f)/2,:)));
            CQT1 = imresize(C_Q_T,[1000,1000]);
            %imagesc(flipud(C_Q_T));
            %figure();
            imagesc(flipud(CQT1));
            imageNameCQT = sprintf("demo_cqt(%d).png",q);
            fullFilePathcqt = fullfile(folderPath, imageNameCQT);
            saveas(gca,fullFilePathcqt)
        end

XTest = [];
Main_Folder ="D:\VK_HMC_Spectrogram_Augmentation_Work\VK_HMC_AUGMENT_FEB_2025\demo\unknown";
count=1;
files1 = dir(Main_Folder);

    for p = 3:length(files1)
        filename = files1(p).name;        
        filename1 = strcat("demo\unknown",'\',filename)
        img = imread(filename1);
        XTest(:,:,1,count) = img(:,:,3);       
        count=count +1;
    end

YTest = classify(augmented_10,XTest);

% Display results
for i = 1:40
    fullFilename = files1(i).name; % Extract filename    
    predictedClass = YTest(i);    
    fprintf('File: %s, Predicted Class: %s\n', fullFilename, string(predictedClass));
end

uniqueStrings = unique(YTest);
counts = zeros(1, length(uniqueStrings));

% Count occurrences of each unique string
for i = 1:length(YTest)
     idx = find(uniqueStrings == YTest(i));
     counts(idx) = counts(idx) + 1;
end

% Find the string with the highest count
[~, maxIndex] = max(counts);
mostFrequentString = uniqueStrings(maxIndex);
fprintf('Prediction by majority voting: %s',mostFrequentString);

toc