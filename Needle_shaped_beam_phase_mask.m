%% Parameter settings
clc;
close all;
clear;

lambda = 1036e-3;
n = 1;
f0 = 300e3;
f = n * f0;

D_OBJ = 17.6;
D_SLM = D_OBJ / 2;
r = D_SLM / 2 * 1e3;

numFoci = 25;
PA = 0.11;
maxIterations = 2;

RL = (lambda / n) * f^2 / (pi * r^2);
focusSpacing = 0.75 * RL;
beamLength = (numFoci - 1) * focusSpacing;
f1 = f - beamLength / 2;
fM = f + beamLength / 2;

dxh = 9.2;
N = round(D_SLM * 1e3 / dxh);
if mod(N, 2) ~= 0
    N = N + 1;
end

slmHeight = 1152;
slmWidth = 1920;
if N > min(slmHeight, slmWidth)
    error('The active mask size N=%d exceeds the SLM dimensions.', N);
end

[x, y] = meshgrid(linspace(-N / 2, N / 2, N) * dxh);
PObj = @(xCoord, yCoord, focalLength) ...
    2 * pi * n .* (sqrt(focalLength.^2 - ...
    (xCoord.^2 + yCoord.^2)) - focalLength) / lambda;

%% Needle-shaped beam phase
pixelAllocation = allocate_pixels(N, numFoci);
gaussianBeam = exp(-(x.^2 + y.^2) / r^2);

[phaseMask, focusPositions, centerIntensity] = optimize_focus_positions( ...
    x, y, PObj, pixelAllocation, gaussianBeam, lambda, n, f0, f, ...
    f1, fM, focusSpacing, PA, maxIterations);

%% Add a linear phase
shiftPixels = 72.5;
kx = (2 * pi / lambda) * (shiftPixels * dxh / f);
ky = kx;
phaseMask = mod(phaseMask + kx * x + ky * y, 2 * pi);

%% Generate, display, and save the SLM phase pattern
phaseArray = uint8(round(phaseMask / (2 * pi) * 255));
PHASE = uint8(zeros(slmHeight, slmWidth) + 128);

startRow = floor((slmHeight - N) / 2) + 1;
endRow = startRow + N - 1;
startCol = floor((slmWidth - N) / 2) + 1;
endCol = startCol + N - 1;
PHASE(startRow:endRow, startCol:endCol) = phaseArray;

figure;
imagesc(PHASE);
axis image;
colormap gray;
colorbar;
title('SLM phase pattern');

folderPath = fullfile(pwd, 'output');
if ~isfolder(folderPath)
    mkdir(folderPath);
end

fileName = sprintf('NB_M%d_%.2fRL_PA%.2f_L%.2gmm_dx%.1f.bmp', ...
    numFoci, focusSpacing / RL, PA, beamLength * 1e-3, shiftPixels);
outputPath = fullfile(folderPath, fileName);
imwrite(PHASE, outputPath);

relativeIntensity = centerIntensity(end, :) / mean(centerIntensity(end, :));
resultsTable = table((1:numFoci).', focusPositions(end, :).', ...
    relativeIntensity.', 'VariableNames', ...
    {'FocusIndex', 'AxialPosition_um', 'RelativeCenterIntensity'});
disp(resultsTable);
fprintf('Phase mask saved to: %s\n', outputPath);


%% Iterative algorithm
function [phaseMask, focusPositions, centerIntensity] = ...
        optimize_focus_positions(x, y, PObj, pixelAllocation, ...
        gaussianBeam, lambda, n, f0, f, f1, fM, focusSpacing, PA, ...
        maxIterations)

    [N, ~, numFoci] = size(pixelAllocation);
    focusPositions = zeros(maxIterations + 1, numFoci);
    centerIntensity = zeros(maxIterations + 1, numFoci);
    focusPositions(1, :) = linspace(f1, fM, numFoci);

    phaseMask = build_phase_mask( ...
        x, y, PObj, pixelAllocation, focusPositions(1, :), f, PA);
    centerIntensity(1, :) = calculate_center_intensity( ...
        phaseMask, focusPositions(1, :), x, y, gaussianBeam, ...
        lambda, n, f0, f, N);

    for iteration = 1:maxIterations
        row = iteration + 1;
        intensityRatio = centerIntensity(row - 1, :) ...
            / mean(centerIntensity(row - 1, :));

        focusPositions(row, 1) = f1;
        for m = 2:numFoci
            if iteration == 1
                previousInterval = focusSpacing;
            else
                previousInterval = focusPositions(row - 1, m) ...
                    - focusPositions(row - 1, m - 1);
            end

            focusPositions(row, m) = focusPositions(row, m - 1) ...
                + previousInterval * intensityRatio(m);
            focusPositions(row, m) = min(max(focusPositions(row, m), f1), fM);
        end
        focusPositions(row, end) = fM;

        phaseMask = build_phase_mask( ...
            x, y, PObj, pixelAllocation, focusPositions(row, :), f, PA);
        centerIntensity(row, :) = calculate_center_intensity( ...
            phaseMask, focusPositions(row, :), x, y, gaussianBeam, ...
            lambda, n, f0, f, N);
    end
end


function phaseMask = build_phase_mask( ...
        x, y, PObj, pixelAllocation, focusPositions, f, PA)

    [N, ~, numFoci] = size(pixelAllocation);
    phaseMask = zeros(N, N);

    for m = 1:numFoci
        focusPhase = PObj(x, y, focusPositions(m)) - PObj(x, y, f);
        phaseMask = phaseMask ...
            + real(focusPhase - PA * m) .* pixelAllocation(:, :, m);
    end

    phaseMask = mod(phaseMask, 2 * pi);
end


function centerIntensity = calculate_center_intensity( ...
        phaseMask, focusPositions, x, y, gaussianBeam, ...
        lambda, n, f0, f, N)

    numFoci = numel(focusPositions);
    centerIntensity = zeros(1, numFoci);
    inputField = gaussianBeam .* exp(1i * phaseMask);
    centerIndex = ceil(N / 2);

    for m = 1:numFoci
        z = abs(focusPositions(m) - f);
        transferFunction = exp(1i * 2 * pi / lambda * z .* ...
            sqrt(1 - (x.^2 + y.^2) * (lambda / n / f0)^2));
        outputField = ifftshift(ifft2(ifftshift( ...
            inputField .* transferFunction)));
        centerIntensity(m) = abs(outputField(centerIndex, centerIndex)).^2;
    end
end


%% Pixel allocation function
function L = allocate_pixels(N, M)

    cellSize = sqrt(M);
    if floor(cellSize) ~= cellSize
        error('The number of foci M must be a perfect square.');
    end

    cellSize = round(cellSize);
    L = zeros(N, N, M);

    for row = 1:cellSize:N
        for col = 1:cellSize:N
            focusOrder = randperm(M);
            index = 1;

            for rowOffset = 0:cellSize - 1
                for colOffset = 0:cellSize - 1
                    if row + rowOffset <= N && col + colOffset <= N
                        L(row + rowOffset, col + colOffset, ...
                            focusOrder(index)) = 1;
                        index = index + 1;
                    end
                end
            end
        end
    end
end
