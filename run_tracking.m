% Main script for SMART tracking. Loads data specified by dataName below,
% and performs SMART tracking, optionally saving a video of the displayed
% output.
% 
% Data can be downloaded using the download_data.m script.
%
% Frank Zijlstra, 2019

clear all
close all

% Add utility functions and NUFFT to path
addpath('./utils');
addpath('./utils/nufft');
addpath('./utils/nufft/utilities');


%% Display options
dispFig = true; % Turn on display of reconstruction
upscaleFactor = 4; % Upscale factor for display
saveVideo = false; % Save display to video file (only works with dispFig = true)


%% Load raw k-space data and template
% dataName = 'spheres_linear';
% % dataName = 'spheres_rotate';
% % dataName = 'spheres_chaos';
% nObjects = 5;

dataName = 'needle_straight';
% dataName = 'needle_45degrees';
nObjects = 1;

load(['./data/data_' dataName])
% Data contains:
% - kLines: [readout X num k-lines X echoes X coils] k-space lines,
% ordered by acquisition time
% - frameIndex: Tells in which fully sampled frame each k-space line was
% acquired
% - profileIndex: Tells which radial profile each k-space line is
% - bitrevFactor: Specifies the bit-reversed sampling factor that was used
% in the acquisition (always 16 in the data provided). This is also the
% acceleration factor
% - frameTime: Acquisition time for each fully sampled frame

load(['./data/template_' dataName])
% Template contains template: [readout X num radial profiles X echoes]


%% Tracking options
dT = frameTime/bitrevFactor;

% System matrix
A = [1 0 dT 0;
     0 1 0 dT;
     0 0 1 0;
     0 0 0 1];
 
% Process noise covariance
Q = diag([dT*1 dT*1 dT*10 dT*10]);

% Observation matrix
H = [1 0 0 0;
     0 1 0 0];
 
% Observation noise covariance
R = 0.5*eye(2);


%%
fprintf('Preparation...\n');

N = max(double(profileIndex)); % Assume maximum profile index is also the image size (true for radial scans on Philips scanners)
nProfiles = N;
nUndersampledProfiles = nProfiles/bitrevFactor;
nReadout = size(kLines,1)/2; % Division by 2 because readout oversampling will be removed

% Calculate radial sampling coordinates
coords = calculateRadialSamplingCoordinates(nReadout, nProfiles);

% Calculate nufft structures for undersampled radial k-space using the
% profileIndex's provided (assumes the pattern repeats bitrevFactor times)
nufftStructure = {};
samplingDensity = {};
for I=1:bitrevFactor
    undersampledCoords = coords(:, profileIndex((I-1)*(nProfiles/bitrevFactor)+1:(I-1)*(nProfiles/bitrevFactor) + nProfiles/bitrevFactor));
   
    ks = [real(undersampledCoords(:)) imag(undersampledCoords(:))]*2*pi;
    nufftStructure{I} = nufft_init(ks, [N N], [6 6], [N N]*2, fftCenter([N N]));
    samplingDensity{I} = bydderSamplingDensity(nufftStructure{I});
end

% Calculate sampling density compensation for the fully sampled radial
% k-space
ks = [real(coords(:)) imag(coords(:))]*2*pi;
samplingDensityFullySampled = bydderSamplingDensity(nufft_init(ks, [N N], [6 6], [N N]*2, fftCenter([N N])));
samplingDensityFullySampled = reshape(samplingDensityFullySampled,nReadout,nProfiles);

%%

% Initialize display
if (dispFig)
    figure(1)
    colormap(gray(256))
    
    set(gcf,'Renderer','zbuffer');
    set(gcf,'Position',[100 100 N*upscaleFactor N*upscaleFactor])
    set(gca,'units','normalized','Position',[0 0 1 1]);
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    hImage = imagesc(zeros(N,N));
    hold on
    plotX = plot(zeros(1,nObjects), zeros(1,nObjects), 'r.', 'MarkerSize',20);
    
    if (saveVideo)
        v = VideoWriter(['movie_' dataName '.avi']); %, 'MPEG-4');
        v.FrameRate = 1 / (frameTime/bitrevFactor);
        v.Quality = 100;
        open(v)
    end
end

% Kalman filter states per object
x = {}; % State vector (x,y,dx,dy)
P = {}; % Uncertainty (4x4 covariance matrix)

initialized = false;
initIm = zeros(N,N);

undersampledKspace = zeros(nReadout, nUndersampledProfiles, size(kLines,3), size(kLines,4));

nufftIndex = 1;
rcsMem = zeros(N,N,size(kLines,4),bitrevFactor);

startTime = tic;

% Loop over all raw k-space profiles
% Note that instead of reading from the stored data, this loop could be
% receiving k-lines from the scanner in real-time.
for I=1:size(kLines,2)
    % Get k-space profile and remove readout oversampling
    undersampledKspace(:,mod(I-1, nUndersampledProfiles)+1,:,:) = removeOversampling(kLines(:,I,:,:),2,1);
    
    % Check if we received a new undersampled frame
    if mod(I-1, nUndersampledProfiles) == nUndersampledProfiles-1
        tic

        % Reconstruct sliding window anatomical image
        rcs = zeros(N,N,size(undersampledKspace,4));
        for C=1:size(undersampledKspace,4)
            rcs(:,:,C) = nufft_adj(double(reshape(undersampledKspace(:,:,1,C) .* samplingDensityFullySampled(:, profileIndex(I-nUndersampledProfiles+1:I)),[],1)), nufftStructure{nufftIndex});
        end
        
        % Store undersampled frame in memory and reconstruct all frames in
        % memory using sum-of-squares coil combination.
        rcsMem(:,:,:,nufftIndex) = rcs;
        anatomicalImage = sqrt(sum(abs(sum(rcsMem,4)).^2,3));
        
        % Reconstruct PC images for all echoes (including the first)
        pcImage = ones(N,N);

        for E=1:size(kLines,3)
            k = undersampledKspace(:,:,E,:);
            
            % Get undersampled template
            undersampledTemplate = template(:,profileIndex(I-nUndersampledProfiles+1:I),E);

            % Flip readouts for even echo numbers
            if (mod(E,2) == 0)
                k = flip(k,1);
            end

            % Perform phase correlation per coil
            rcs = zeros(N,N,size(k,4));
            for C=1:size(k,4)
                % Perform phase correlation
                M = k(:,:,:,C) .* conj(undersampledTemplate);
                M = sign(M);

                rcs(:,:,C) = nufft_adj(double(M(:)) .* samplingDensity{nufftIndex}(:), nufftStructure{nufftIndex});
            end

            % Sum of squares reconstruction and multiplication with
            % combined PC image
            pcImage = pcImage .* abs(sqrt(sum(abs(rcs).^2,3)));
        end
        
        % Increment nufftStructure index 
        nufftIndex = mod(nufftIndex,bitrevFactor)+1;
        
        
        % Display sliding window anatomical image
        if (dispFig)
            set(hImage, 'CData', abs(anatomicalImage));
        end

        
        % Perform tracking...

        % Wait for first full frame to be acquired
        if (frameIndex(I) == 1)
            % Combine undersampled PC images in first fully sampled image
            % for initialization. In the ideal case the reconstruction
            % would be performed with different weights for the fully
            % sampled recon. But this appears to work.
            initIm = initIm + abs(pcImage);
            continue;
        end

        % If first full frame is acquired, initialize object positions
        if (~initialized && frameIndex(I) == 2)
            for O=1:nObjects
                % Find highest intensity match coordinates
                [~,ind] = max(initIm(:));
                [bestMatchY,bestMatchX] = ind2sub(size(initIm), ind);

                % Mask out neighbouring matches
                RD = 7;
                initIm(bestMatchY-RD:bestMatchY+RD,bestMatchX-RD:bestMatchX+RD) = 0;

                % Enter coordinates into Kalman filter's initial state
                x{O} = [bestMatchX bestMatchY 0 0].';

                % Initialize Kalman filter's uncertainty
                P{O} = 2*eye(4);
            end
            
            initialized = true;
            continue;
        end

        % Kalman filter prediction step
        for O=1:nObjects
            x{O} = A*x{O};
            P{O} = A*P{O}*A' + Q;
        end

        % Search for matches in PC image
        searchImage = abs(pcImage);
        
        % Determine noise level
        stdIm = std(searchImage(:));
        
        % Keep only local maxima and filter matches below noise level
        searchImage = (searchImage==imdilate(searchImage,ones(7,7))).*searchImage;
        searchImage(searchImage < 1*stdIm) = 0;

        % Keep only 50 coordinates with highest PC correlation
        [~,matchIndices] = sort(searchImage(:),'descend');
        matchIndices = matchIndices(1:min(size(matchIndices,1),50));
        
        [my, mx] = ind2sub(size(searchImage),matchIndices);
        pcCorrelations = searchImage(matchIndices);

        % Rank matches by distance to current position divided by PC value
        % (i.e. find a close, high intensity match)
        currentPositions = cell2mat(x);
        currentPositions = currentPositions(1:2,:); % Keep only x and y coordinate
        D = (pdist2(currentPositions.',[mx my])+1).^1 ./ repmat((pcCorrelations.'),size(currentPositions,2),1);
        matchIndices = munkres(D);
        
        % Kalman filter update step
        for O=1:nObjects
            z = [mx(matchIndices(O)) my(matchIndices(O))].';

            K = P{O}*H'*inv(H*P{O}*H' + R);
            x{O} = x{O} + K*(z - H*x{O});
            P{O} = P{O} - K*H*P{O};
        end
        
        % Display new object locations
        if (dispFig)
            tmp = cell2mat(x);
            set(plotX, 'XData', tmp(1,:), 'YData', tmp(2,:));
            drawnow
            
            if (saveVideo)
                VF = getframe;
                writeVideo(v, VF.cdata)
            end
        end

        dur = toc;
        fprintf('FPS: %.1f\n', 1/dur);
    end
end

fprintf('Total time: %.1f seconds (%.1f FPS)\n', toc(startTime), (size(kLines,2)/nUndersampledProfiles) / toc(startTime));

%% Close video
if (dispFig)
    if (saveVideo)
        close(v);
    end
end
