
function startSimulationSen()
clc;
clear all;
close all;
% MOANS DV Hop andAS-IDV-Hop
numAnchorsSen = 15;
numUnknownNodesSen = 85;
communicationRangeSen = 30;
maxOptimizedAnchors = 15; 
 maxIterationsLoop = 6; 

anchorPositionsSen = rand(numAnchorsSen, 2) * 100;
unknownNodePositionsSen = rand(numUnknownNodesSen, 2) * 100;

figure;
hold on;
scatter(anchorPositionsSen(:, 1), anchorPositionsSen(:, 2), 100, 'r*');
scatter(unknownNodePositionsSen(:, 1), unknownNodePositionsSen(:, 2), 100, 'bo');
title('Node Positions');
xlabel('X-axis');
ylabel('Y-axis');
legend('All Anchors', 'Unknown Nodes');
grid on;
axis equal;
numSelect=maxOptimizedAnchors;
Num_Solutions=50;
Max_Iter=100;
unknownNodePositions=unknownNodePositionsSen;
anchorPositions=anchorPositionsSen;

optimizedAnchorIndices = selectAnchorNodesEO(anchorPositions, numSelect, unknownNodePositions, Max_Iter, Num_Solutions);

optimizedAnchorPositions = anchorPositionsSen(optimizedAnchorIndices, :);


scatter(optimizedAnchorPositions(:, 1), optimizedAnchorPositions(:, 2), 200, 'filled', 'MarkerFaceColor', 'g');
legend('All Anchors', 'Unknown Nodes', 'Selected Anchors');
hold off;


Max_Iter = 100;
Number_of_particles = 50;
Threshold = 0.5;

best_candidates = zeros(4, numUnknownNodesSen, 2);
pbest = Inf;

x_i = generatePositionVector(optimizedAnchorPositions);


for i = 1:Max_Iter
    for j = 1:Number_of_particles
        estimated_coords = LSMestimate(x_i, optimizedAnchorPositions, communicationRangeSen, numUnknownNodesSen);
        fitness = evaluateFitness(estimated_coords, unknownNodePositionsSen);
        [best_candidates, pbest] = updateBestCandidates(fitness, estimated_coords, best_candidates, pbest);
        X_eq = computeEquilibrium(best_candidates);
        numParticles = size(x_i, 1);
        dim = size(x_i, 2);
        F = calculateF(X_eq, numParticles, dim);
        G = calculateG(X_eq, numParticles, dim);
        x_i = updatePosition(x_i, F, G);
    end
end

OANS = best_candidates(1, :, :);
X_binary = (OANS > Threshold);


hop_size = recalculateHopSize(OANS, communicationRangeSen);





function indices = selectAnchorNodesEO(anchorPositions, numSelect, unknownNodePositions, Max_Iter, Num_Solutions)


    solutions = rand(Num_Solutions, size(anchorPositions, 1)) > 0.5;
    
   
    for i = 1:Num_Solutions
        currentSolution = find(solutions(i, :));
        if numel(currentSolution) > numSelect
            currentSolution = currentSolution(1:numSelect);
        end
        solutions(i, :) = false;
        solutions(i, currentSolution) = true;
    end

    bestFitness = Inf;
    bestIndices = [];

 
    for iter = 1:Max_Iter
        for i = 1:Num_Solutions
            indices = find(solutions(i, :));
            fitness = evaluateAnchorFitness(anchorPositions(indices, :), unknownNodePositions);
            
           
            if fitness < bestFitness
                bestFitness = fitness;
                bestIndices = indices;
            end
        end
        
        
        X_eq = mean(solutions, 1);
        for i = 1:Num_Solutions
            newSol = rand(1, size(anchorPositions, 1)) < X_eq;
            newIndices = find(newSol);
            if numel(newIndices) > numSelect
                newIndices = newIndices(1:numSelect);
            end
            solutions(i, :) = false;
            solutions(i, newIndices) = true;
        end
    end

    indices = bestIndices;
end

function fitness = evaluateAnchorFitness(selectedAnchorPositions, unknownNodePositions)
   
    totalDistance = 0;
    for i = 1:size(unknownNodePositions, 1)
        distances = vecnorm(unknownNodePositions(i, :) - selectedAnchorPositions, 2, 2);
        totalDistance = totalDistance + min(distances);
    end
    fitness = totalDistance / size(unknownNodePositions, 1);
end


function x_i = generatePositionVector(anchorPositions)
  
    numParticles = 50;
    numAnchors = size(anchorPositions, 1);
    dim = size(anchorPositions, 2); 
    expandedAnchorPositions = repmat(anchorPositions, numParticles, 1);
    randomDisplacements = (rand(numParticles * numAnchors, dim) - 0.5) * 0.1;

    x_i = expandedAnchorPositions + randomDisplacements;
end

function estimated_coords = LSMestimate(particlePositions, anchorPositions, communicationRange, numUnknownNodes)
   
    estimated_coords = zeros(numUnknownNodes, 2);

    for unknownIdx = 1:numUnknownNodes
        distances = vecnorm(anchorPositions - repmat(particlePositions(unknownIdx, :), size(anchorPositions, 1), 1), 2, 2);
        inRangeIndices = find(distances <= communicationRange);
        if length(inRangeIndices) < 3
            inRangeIndices = 1:size(anchorPositions, 1);
        end
        anchorsInRange = anchorPositions(inRangeIndices, :);
        distancesInRange = distances(inRangeIndices);
        estimated_coords(unknownIdx, :) = mean(anchorsInRange, 1);
    end
end

function fitness = evaluateFitness(estimated_coords, unknownNodePositions)
    differences = unknownNodePositions - estimated_coords;
    fitness = sum(vecnorm(differences, 2, 2)); 
end

function [best_candidates, pbest] = updateBestCandidates(fitness, current_coords, best_candidates, pbest)
    if fitness < pbest
        pbest = fitness;
        best_candidates(1, :, :) = current_coords;
    end
end

function X_eq = computeEquilibrium(best_candidates)
    X_eq = mean(best_candidates, 1);  
    X_eq = squeeze(X_eq); 
end

function F = calculateF(X_eq, numParticles, dim)
    F = rand(numParticles, dim);
end

function G = calculateG(X_eq, numParticles, dim)
    G = rand(numParticles, dim);
end

function new_position = updatePosition(current_position, F, G)
    new_position = current_position + F - G;
end

function hop_size = recalculateHopSize(OANS, communicationRangeSen)
    hop_size = mean(OANS) * communicationRangeSen;
end

numAnchorsSen=size(optimizedAnchorPositions, 1);
numUnknownNodesSen=numUnknownNodesSen;
anchorPositionsSen=optimizedAnchorPositions;
unknownNodePositionsSen=unknownNodePositionsSen;
communicationRangeSen=communicationRangeSen;
disp(['numAnchorsSen: ', num2str(numAnchorsSen)]);
disp(['numUnknownNodesSen: ', num2str(numUnknownNodesSen)]);
RSSISen = generateRSSIValuesSen(numAnchorsSen, numUnknownNodesSen, anchorPositionsSen, unknownNodePositionsSen,communicationRangeSen); 
     minHopCountSen = findMinHopCountSen(RSSISen, numAnchorsSen, numUnknownNodesSen);
     
     averageHopDistanceSen = calculateAverageHopDistance(anchorPositionsSen, unknownNodePositionsSen, minHopCountSen);
     estimatedUnknownNodePositionsSen = calculateUnknownNodeCoordinates(anchorPositionsSen, minHopCountSen,averageHopDistanceSen);

    errors = sqrt(sum(((unknownNodePositionsSen - estimatedUnknownNodePositionsSen)).^2, 2));
    figure;
   hold on;
    plot(1:numUnknownNodesSen, errors, '-og');
    title('Accuracy Graph for Each Unknown Node');
    xlabel('Unknown Node Index');
    ylabel('Positioning Error');
    grid on;

sizeUpdated = size(unknownNodePositionsSen);
sizeEstimated = size(estimatedUnknownNodePositionsSen);
if isequal(sizeUpdated, sizeEstimated)
    differences_dv = abs(unknownNodePositionsSen - estimatedUnknownNodePositionsSen);
    thresholdSen = 0.9;
    accuratePredictions = differences_dv < thresholdSen;
    numAccuratePredictions = sum(accuratePredictions(:));
    accuracy_dv = numAccuratePredictions / numel(unknownNodePositionsSen);
    accuracy_val_dv=num2str(accuracy_dv * 100);
    disp(accuracy_val_dv);
    disp(['Accuracy: ', num2str(accuracy_dv * 100), '%']);
else
    disp('Error: Matrix dimensions must agree.');
end

    function RSSISen = generateRSSIValuesSen(numAnchorsSen, numUnknownNodesSen, anchorPositionsSen, unknownNodePositionsSen,communicationRangeSen)
    P_trSen = 20; 
    P_lossSen = 30;
    d_0Sen = 1; 
    tauSen = 2; 
    RSSISen = zeros(numAnchorsSen, numUnknownNodesSen);

    for i = 1:numAnchorsSen
        for j = 1:numUnknownNodesSen
            distance = sqrt((anchorPositionsSen(i, 1) - unknownNodePositionsSen(j, 1))^2 + (anchorPositionsSen(i, 2) - unknownNodePositionsSen(j, 2))^2);
            RSSISen(i, j) = P_trSen - P_lossSen - 10 * tauSen * log10(distance / d_0Sen);
        end
    end
end

function minHopCountSen = findMinHopCountSen(RSSISen, numAnchorsSen, numUnknownNodesSen)
    minHopCountSen = zeros(numAnchorsSen, numUnknownNodesSen);
    for j = 1:numUnknownNodesSen
        hopCountSen = inf(numAnchorsSen, 1); 

        [~, anchorIdxSen] = max(RSSISen(:, j));
        hopCountSen(anchorIdxSen) = 0; 


        for i = 1:numAnchorsSen
            for k = 1:numAnchorsSen
                if RSSISen(k, j) > RSSISen(i, j)
                    hopCountSen(i) = min(hopCountSen(i), hopCountSen(k) + 1);
                end
            end
        end

    
        minHopCountSen(:, j) = hopCountSen;
    end
end

  function averageHopDistanceSen = calculateAverageHopDistance(anchorPositionsSen, unknownNodePositionsSen, minHopCountSen)
    numAnchorsSen = size(anchorPositionsSen, 1);
    numUnknownNodesSen = size(unknownNodePositionsSen, 1);
    
    if size(minHopCountSen, 1) < numAnchorsSen || size(minHopCountSen, 2) < numUnknownNodesSen
        error('minHopCountSen does not have the correct dimensions.');
    end

    averageHopDistanceSen = zeros(numAnchorsSen, 1);

    for i = 1:numAnchorsSen
        totalCorrectedDistanceSen = 0;
        totalCountSen = 0;

        for j = 1:numUnknownNodesSen
            if minHopCountSen(i, j) > 0
                actualDistance = sqrt((anchorPositionsSen(i, 1) - unknownNodePositionsSen(j, 1))^2 + (anchorPositionsSen(i, 2) - unknownNodePositionsSen(j, 2))^2);
                correctedDistance = actualDistance / minHopCountSen(i, j);
                totalCorrectedDistanceSen = totalCorrectedDistanceSen + correctedDistance;
                totalCountSen = totalCountSen + 1;
            end
        end

        if totalCountSen > 0
            averageHopDistanceSen(i) = totalCorrectedDistanceSen / totalCountSen;
        end
    end
end
function estimatedUnknownNodePositionsSen = calculateUnknownNodeCoordinates(anchorPositionsSen, minHopCountSen, averageHopDistanceSen)
    numUnknownNodesSen = size(minHopCountSen, 2);
    estimatedUnknownNodePositionsSen = zeros(numUnknownNodesSen, 2);

    for j = 1:numUnknownNodesSen
        sumXSen = 0;
        sumYSen = 0;
        totalHopDistance = 0;

        for i = 1:size(anchorPositionsSen, 1)
            sumXSen = sumXSen + anchorPositionsSen(i, 1) * minHopCountSen(i, j) / averageHopDistanceSen(i);
            sumYSen = sumYSen + anchorPositionsSen(i, 2) * minHopCountSen(i, j) / averageHopDistanceSen(i);
            totalHopDistance = totalHopDistance + minHopCountSen(i, j) / averageHopDistanceSen(i);
        end

        if totalHopDistance > 0
            estimatedUnknownNodePositionsSen(j, 1) = sumXSen / totalHopDistance;
            estimatedUnknownNodePositionsSen(j, 2) = sumYSen / totalHopDistance;
        end
    end
end

 numAnchors_updated_itrSen=numAnchorsSen;
    numUnknownNodes_updated_itrSen=numUnknownNodesSen;
    anchorPositions_updated_itrSen=anchorPositionsSen;
    unknownNodePositions_updated_itrSen=unknownNodePositionsSen;
    

    disp(anchorPositions_updated_itrSen);
    updatedAnchorPositions = removeCollinearAnchors(anchorPositions_updated_itrSen);

disp(updatedAnchorPositions);
anchorPositions_updated_itrSen=updatedAnchorPositions;
numAnchors_updated_itrSen=length(updatedAnchorPositions);
function updatedAnchorPositionsSen = removeCollinearAnchors(anchorPositions_updated_itrSen)
    updatedAnchorPositionsSen = anchorPositions_updated_itrSen;
    i = 1;

    while i <= size(updatedAnchorPositionsSen, 1) - 2
        A = updatedAnchorPositionsSen(i, :);
        for j = i + 1:size(updatedAnchorPositionsSen, 1) - 1
            B = updatedAnchorPositionsSen(j, :);
            for k = j + 1:size(updatedAnchorPositionsSen, 1)
                C = updatedAnchorPositionsSen(k, :);
                if isCollinearTripletSen(A, B, C)
                    updatedAnchorPositionsSen(k, :) = [];
                    break; 
                end
            end
            if size(updatedAnchorPositionsSen, 1) < size(anchorPositionsSen, 1)
                break; 
            end
        end
        if size(updatedAnchorPositionsSen, 1) < size(anchorPositionsSen, 1)
            continue; 
        end
        i = i + 1;
    end
end

function isCollinear = isCollinearTripletSen(A, B, C)
    area = abs(0.5 * (A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2))));
    isCollinear = area < 1e-6; 
end
    
    
imp_accuracySen = cell(maxIterationsLoop, 1);
for iterationSen = 1:maxIterationsLoop 
   disp('iterations');
   disp(iterationSen); 
 RSSI_itr = generateRSSIValues_itr(numAnchors_updated_itrSen, numUnknownNodes_updated_itrSen, anchorPositions_updated_itrSen, unknownNodePositions_updated_itrSen,communicationRangeSen);

 % figure;
  minHopCount_itr = findMinHopCount_itr(RSSI_itr, numAnchors_updated_itrSen, numUnknownNodes_updated_itrSen);
   disp(minHopCount_itr);
    IMCDS_itr = ImprovedMinimumConnectedDominationSet_itr(numAnchors_updated_itrSen, numUnknownNodes_updated_itrSen, anchorPositions_updated_itrSen, unknownNodePositions_updated_itrSen, minHopCount_itr, RSSI_itr);
    disp(['iterationSen ', num2str(iterationSen), ' IMCDS_itr:']);
    disp('Improved Minimum Connected Dominating Set (Anchor Indices):');
    
     averageHopDistance_again_itr = calculateAverageHopDistance_again_itr(IMCDS_itr, unknownNodePositions_updated_itrSen, minHopCount_itr); 
    estimatedUnknownNodePositions_again_itr = calculateUnknownNodeCoordinates_again_itr(IMCDS_itr, minHopCount_itr, averageHopDistance_again_itr)
  
       figure;
    errors_again_itr = sqrt(sum(((unknownNodePositions_updated_itrSen - estimatedUnknownNodePositions_again_itr)).^2, 2));
    disp('Positional Errors for Each Unknown Node again itr:');
    disp(errors_again_itr);
    hold on;
color55 = [0.4, 0.4, 0.4];

      plot(1:numUnknownNodesSen, errors, '-o', 'Color', color55);
     plot(1:numUnknownNodes_updated_itrSen, errors_again_itr, '-ob');
    title('Accuracy Graph for Each Unknown Node');
    xlabel('Unknown Node Index');
    ylabel('Positioning Error');
    grid on;
   
sizeUpdated = size(unknownNodePositions_updated_itrSen);
sizeEstimated = size(estimatedUnknownNodePositions_again_itr);
if isequal(sizeUpdated, sizeEstimated)
    differences = abs(unknownNodePositions_updated_itrSen - estimatedUnknownNodePositions_again_itr);
    thresholdSen = 0.9;
   
    accuratePredictions = differences < thresholdSen;
    numAccuratePredictions = sum(accuratePredictions(:));
    accuracy = numAccuratePredictions / numel(unknownNodePositions_updated_itrSen);
    accuracy_val=num2str(accuracy * 100);
    disp(accuracy_val);
    
    disp(['Accuracy: ', num2str(accuracy * 100), '%']);
else
    disp('Error: Matrix dimensionss.');
end

  imp_accuracySen{iterationSen} = accuracy_val;
 
    nodeIndices_itrSen = (1:size(errors_again_itr, 1))';
    extendedData_itrSen = [nodeIndices_itrSen, errors_again_itr];
    sortedData_itr = sortrows(extendedData_itrSen, 2);       
    lowestLocalizationError_itr = sortedData_itr(1:5, 1:2);
   lowest_error_list_itr = lowestLocalizationError_itr(:, 1);
  
if exist('unknownNodePositions_updated_itrSen', 'var')
    validIndices = find(lowest_error_list_itr <= size(unknownNodePositions_updated_itrSen, 1));
    if ~isempty(validIndices)
        removedRows = unknownNodePositions_updated_itrSen(validIndices, :);
        unknownNodePositions_updated_itrSen(validIndices, :) = [];
       
        
    else
        disp('No rows to remove.');
    end
else

    disp('Error: unknownNodePositions_updated_itrSen is not defined.');

    break; 
end
    assistantAnchorNodes_itrSen = [];
    assistantAnchorNodes_itrSen = [assistantAnchorNodes_itrSen; removedRows];
    numAnchors_updated_itrSen = numAnchors_updated_itrSen + length(assistantAnchorNodes_itrSen);
   
   numUnknownNodes_updated_itrSen=length(unknownNodePositions_updated_itrSen);
     anchorPositions_updated_itrSen = [anchorPositions_updated_itrSen; assistantAnchorNodes_itrSen];

   
   
end




numericArray = cellfun(@str2double, imp_accuracySen);


disp(numericArray);

disp(accuracy_val_dv);




end



function RSSI_itr = generateRSSIValues_itr(numAnchorsSen, numUnknownNodesSen, anchorPositionsSen, unknownNodePositionsSen,communicationRangeSen)
disp('numAnchorsSen');
disp(numAnchorsSen);
disp('numUnknownNodesSen');
disp(numUnknownNodesSen);
disp(anchorPositionsSen);
    P_trSen = 20; 
    P_lossSen = 30; 
    d_0Sen = 1; 
    tauSen = 2; 
    RSSI_itr = zeros(numAnchorsSen, numUnknownNodesSen);

    for i = 1:numAnchorsSen
        for j = 1:numUnknownNodesSen
            distance = sqrt((anchorPositionsSen(i, 1) - unknownNodePositionsSen(j, 1))^2 + (anchorPositionsSen(i, 2) - unknownNodePositionsSen(j, 2))^2);
            RSSI_itr(i, j) = P_trSen - P_lossSen - 10 * tauSen * log10(distance / d_0Sen);
        end
    end
end
 
function minHopCount_itr = findMinHopCount_itr(RSSI_itr, numAnchors_updated_itrSen, numUnknownNodes_updated_itrSen)
    minHopCount_itr = zeros(numAnchors_updated_itrSen, numUnknownNodes_updated_itrSen);

    for j = 1:numUnknownNodes_updated_itrSen
        hopCount_itr = inf(numAnchors_updated_itrSen, 1); 

        [~, anchorIdx_itr] = max(RSSI_itr(:, j));
        hopCount_itr(anchorIdx_itr) = 0; 

        for i = 1:numAnchors_updated_itrSen
            for k = 1:numAnchors_updated_itrSen
                if RSSI_itr(k, j) > RSSI_itr(i, j)
                    hopCount_itr(i) = min(hopCount_itr(i), hopCount_itr(k) + 1);
                end
            end
        end

        minHopCount_itr(:, j) = hopCount_itr;
    end
end


function IMCDS_itr = ImprovedMinimumConnectedDominationSet_itr(numAnchorsSen, numUnknownNodesSen, anchorPositionsSen, unknownNodePositionsSen, minHopCountSen, RSSISen)

    IMCDS_Indices = [];

    minRequiredAnchors = ceil(numAnchorsSen / 3);
    
    connectivity = sum(minHopCountSen ~= inf, 2); 
    [~, sortedIndices] = sort(connectivity, 'descend'); 

    coveredUnknownNodes = false(numUnknownNodesSen, 1);
    for i = 1:numAnchorsSen
        anchorIdxSen = sortedIndices(i);
        connectedUnknowns = minHopCountSen(anchorIdxSen, :) ~= inf;
        isNewNodeCovered = any(~coveredUnknownNodes & connectedUnknowns);
        isBelowMinThreshold = length(IMCDS_Indices) < minRequiredAnchors;
        isNewNodeCovered = isscalar(isNewNodeCovered) && isNewNodeCovered;
        isBelowMinThreshold = isscalar(isBelowMinThreshold) && isBelowMinThreshold;
        if isNewNodeCovered || isBelowMinThreshold
            IMCDS_Indices = [IMCDS_Indices, anchorIdxSen];
            coveredUnknownNodes = coveredUnknownNodes | connectedUnknowns;
        end
allNodesCovered = isscalar(all(coveredUnknownNodes))&& coveredUnknownNodes;
       
            if((allNodesCovered &&length(IMCDS_Indices))>= minRequiredAnchors)
            break;
        end
    end

    IMCDS_itr = anchorPositionsSen(IMCDS_Indices, :);
end

function averageHopDistance_again_itr = calculateAverageHopDistance_again_itr(anchorPositionsSen, unknownNodePositionsSen, minHopCountRecalculate)
    numAnchorsSen = size(anchorPositionsSen, 1);
    disp(minHopCountRecalculate);
    numUnknownNodesSen = size(unknownNodePositionsSen, 1);
    if size(minHopCountRecalculate, 1) < numAnchorsSen || size(minHopCountRecalculate, 2) < numUnknownNodesSen
        error('minHopCountSen does not have the correct dimensions.');
    end

    averageHopDistance_again_itr = zeros(numAnchorsSen, 1);

    for i = 1:numAnchorsSen
        totalCorrectedDistanceSen = 0;
        totalCountSen = 0;

        for j = 1:numUnknownNodesSen
            if minHopCountRecalculate(i, j) > 0
                actualDistance = sqrt((anchorPositionsSen(i, 1) - unknownNodePositionsSen(j, 1))^2 + (anchorPositionsSen(i, 2) - unknownNodePositionsSen(j, 2))^2);
                correctedDistance = actualDistance / minHopCountRecalculate(i, j); 
                totalCorrectedDistanceSen = totalCorrectedDistanceSen + correctedDistance;
                totalCountSen = totalCountSen + 1;
            end
        end

        if totalCountSen > 0
            averageHopDistance_again_itr(i) = totalCorrectedDistanceSen / totalCountSen;
        end
    end
end

function estimatedUnknownNodePositions_again_itr = calculateUnknownNodeCoordinates_again_itr(anchorPositionsSen, minHopCountSen, averageHopDistanceSen)
    numUnknownNodesSen = size(minHopCountSen, 2);
    estimatedUnknownNodePositions_again_itr = zeros(numUnknownNodesSen, 2);

    for j = 1:numUnknownNodesSen
        sumXSen = 0;
        sumYSen = 0;
        tHopDistance = 0;

        for i = 1:size(anchorPositionsSen, 1)
            sumXSen = sumXSen + anchorPositionsSen(i, 1) * minHopCountSen(i, j) / averageHopDistanceSen(i);
            sumYSen = sumYSen + anchorPositionsSen(i, 2) * minHopCountSen(i, j) / averageHopDistanceSen(i);
            tHopDistance = tHopDistance + minHopCountSen(i, j) / averageHopDistanceSen(i);
        end
        if tHopDistance > 0
            estimatedUnknownNodePositions_again_itr(j, 1) = sumXSen / tHopDistance;
            estimatedUnknownNodePositions_again_itr(j, 2) = sumYSen / tHopDistance;
        end
    end
    end