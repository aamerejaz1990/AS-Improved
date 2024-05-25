
function startSimulation()
clc;
% MCDS-DVHOP-ANDAS-IDV-Hop
clear all;
close all;
    disp('Simulation started.');

    
    numAnchors = 20;
    numUnknownNodes = 80;
    communicationRange = 30;
     maxIterations = 6; 

    anchorPositions = rand(numAnchors, 2) * 100;
    unknownNodePositions = rand(numUnknownNodes, 2) * 100;
   
    color1=[31,120,180]/255;
    color2 = [49,163,84] / 255;
    figure;
    hold on;
    scatter(anchorPositions(:, 1), anchorPositions(:, 2), 100, 'r*'); % Anchor nodes
    scatter(unknownNodePositions(:, 1), unknownNodePositions(:, 2), 100, 'bo'); % Unknown nodes
    title('Node Positions');
    xlabel('X-axis');
    ylabel('Y-axis');
    legend('Anchors', 'Unknown Nodes');
    grid on;
    axis equal;
    hold off;


      RSSI = generateRSSIValues(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions,communicationRange);
  
    minHopCount = findMinHopCount(RSSI, numAnchors, numUnknownNodes);
 

    IMCDS_Coordinates = ImprovedMinimumConnectedDominationSet(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions, minHopCount, RSSI);
  
 
    
    numAnchorss = length(IMCDS_Coordinates);

     RSSI_Updated = generateRSSIValuesUpdated(IMCDS_Coordinates, numUnknownNodes, unknownNodePositions,communicationRange);

minHopCountRecalculate = findMinHopCountRecalculate(RSSI_Updated, numAnchorss, numUnknownNodes);

averageHopDistance_again = calculateAverageHopDistance_again(IMCDS_Coordinates, unknownNodePositions, minHopCountRecalculate);

estimatedUnknownNodePositions_again = calculateUnknownNodeCoordinates_again(IMCDS_Coordinates, minHopCountRecalculate,averageHopDistance_again);

    errors_again = sqrt(sum(((unknownNodePositions - estimatedUnknownNodePositions_again)).^2, 2));
  
     figure;
     hold on;
    plot(1:numUnknownNodes, errors_again, '-.o', 'Color', color1, 'LineWidth', 2, 'MarkerSize', 8);
    title('Accuracy Graph for Each Unknown Node');
    xlabel('Unknown Node Index');
    ylabel('Positioning Error');
     legend('MCDS-DV-HOP', 'AS-IDV-Hop');
    grid on;
    
    
        
sizeUpdated = size(unknownNodePositions);
sizeEstimated = size(estimatedUnknownNodePositions_again);

if isequal(sizeUpdated, sizeEstimated)
  
    differences_dv = abs(unknownNodePositions - estimatedUnknownNodePositions_again);
    
    threshold = 0.9;
    accuratePredictions = differences_dv < threshold;
    numAccuratePredictions = sum(accuratePredictions(:));
    accuracy_dv = numAccuratePredictions / numel(unknownNodePositions);
    accuracy_val_dv=num2str(accuracy_dv * 100);
    
else
    
    disp('Error');
end

    
   function RSSI = generateRSSIValues(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions,communicationRange)
    P_tr = 20; 
    P_loss = 30; 
    d_0 = 1; 
    tau = 2; 
    RSSI = zeros(numAnchors, numUnknownNodes);

    for i = 1:numAnchors
        for j = 1:numUnknownNodes
            distance = sqrt((anchorPositions(i, 1) - unknownNodePositions(j, 1))^2 + (anchorPositions(i, 2) - unknownNodePositions(j, 2))^2);
            RSSI(i, j) = P_tr - P_loss - 10 * tau * log10(distance / d_0);
        end
    end
   end



function minHopCount = findMinHopCount(RSSI, numAnchors, numUnknownNodes)

    minHopCount = zeros(numAnchors, numUnknownNodes);

    for j = 1:numUnknownNodes
        hopCount = inf(numAnchors, 1); 

    
        [~, anchorIdx] = max(RSSI(:, j));
        hopCount(anchorIdx) = 0; 

        for i = 1:numAnchors
            for k = 1:numAnchors
                if RSSI(k, j) > RSSI(i, j)
                    hopCount(i) = min(hopCount(i), hopCount(k) + 1);
                end
            end
        end

        minHopCount(:, j) = hopCount;
    end
end


function IMCDS_Coordinates = ImprovedMinimumConnectedDominationSet(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions, minHopCount, RSSI)
    
    IMCDS_Indices = [];

  
    minRequiredAnchors = ceil(numAnchors / 3);

    connectivity = sum(minHopCount ~= inf, 2); 
    [~, sortedIndices] = sort(connectivity, 'descend');

    coveredUnknownNodes = false(numUnknownNodes, 1);
    for i = 1:numAnchors
        anchorIdx = sortedIndices(i);
        connectedUnknowns = minHopCount(anchorIdx, :) ~= inf;

        isNewNodeCovered = any(~coveredUnknownNodes & connectedUnknowns);
        
       
        isBelowMinThreshold = length(IMCDS_Indices) < minRequiredAnchors;

        isNewNodeCovered = isscalar(isNewNodeCovered) && isNewNodeCovered;
        isBelowMinThreshold = isscalar(isBelowMinThreshold) && isBelowMinThreshold;
        if isNewNodeCovered || isBelowMinThreshold
            IMCDS_Indices = [IMCDS_Indices, anchorIdx];
            coveredUnknownNodes = coveredUnknownNodes | connectedUnknowns;
        end
allNodesCovered = isscalar(all(coveredUnknownNodes))&& coveredUnknownNodes;
 
       
       
            if((allNodesCovered &&length(IMCDS_Indices))>= minRequiredAnchors)
            break;
        end
    end

    
    IMCDS_Coordinates = anchorPositions(IMCDS_Indices, :);
end



function RSSI_Updated = generateRSSIValuesUpdated(IMCDS, numUnknownNodes, unknownNodePositions,communicationRange)
    P_tr = 20; 
    P_loss = 30; 
    d_0 = 1; 
    tau = 2; 
    RSSI_Updated = zeros(length(IMCDS), numUnknownNodes);

    for i = 1:length(IMCDS)
        for j = 1:numUnknownNodes
            distance = sqrt((IMCDS(i, 1) - unknownNodePositions(j, 1))^2 + (IMCDS(i, 2) - unknownNodePositions(j, 2))^2);
            RSSI_Updated(i, j) = P_tr - P_loss - 10 * tau * log10(distance / d_0);
        end
    end
end

function minHopCountRecalculate = findMinHopCountRecalculate(RSSI_Updated, numAnchorss, numUnknownNodes)
  
    minHopCountRecalculate = zeros(numAnchorss, numUnknownNodes);
    for j = 1:numUnknownNodes
        hopCountRecal = inf(numAnchorss, 1); 

        [~, anchorIdx] = max(RSSI_Updated(:, j));
        hopCountRecal(anchorIdx) = 0; 
        for i = 1:numAnchorss
            for k = 1:numAnchorss
               if RSSI_Updated(k, j) > RSSI_Updated(i, j)
                    hopCountRecal(i) = min(hopCountRecal(i), hopCountRecal(k) + 1);
                end
            end
       end

        minHopCountRecalculate(:, j) = hopCountRecal;
    end
end


function averageHopDistance_again = calculateAverageHopDistance_again(anchorPositions, unknownNodePositions, minHopCountRecalculate)

    numAnchors = size(anchorPositions, 1);
    disp(minHopCountRecalculate);
    numUnknownNodes = size(unknownNodePositions, 1);

    if size(minHopCountRecalculate, 1) < numAnchors || size(minHopCountRecalculate, 2) < numUnknownNodes
        error('minHopCount does not have the correct dimensions.');
    end

    averageHopDistance_again = zeros(numAnchors, 1);

    for i = 1:numAnchors
        totalCorrectedDistance = 0;
        totalCount = 0;

        for j = 1:numUnknownNodes
            if minHopCountRecalculate(i, j) > 0
                actualDistance = sqrt((anchorPositions(i, 1) - unknownNodePositions(j, 1))^2 + (anchorPositions(i, 2) - unknownNodePositions(j, 2))^2);
                correctedDistance = actualDistance / minHopCountRecalculate(i, j); % Correction factor
                totalCorrectedDistance = totalCorrectedDistance + correctedDistance;
                totalCount = totalCount + 1;
            end
        end

        if totalCount > 0
            averageHopDistance_again(i) = totalCorrectedDistance / totalCount;
        end
    end
end

function estimatedUnknownNodePositions_again = calculateUnknownNodeCoordinates_again(anchorPositions, minHopCount, averageHopDistance)
    numUnknownNodes = size(minHopCount, 2);
    estimatedUnknownNodePositions_again = zeros(numUnknownNodes, 2);

    for j = 1:numUnknownNodes
        % Initialize 
        sumX = 0;
        sumY = 0;
        totalHopDistance = 0;

        for i = 1:size(anchorPositions, 1)
            sumX = sumX + anchorPositions(i, 1) * minHopCount(i, j) / averageHopDistance(i);
            sumY = sumY + anchorPositions(i, 2) * minHopCount(i, j) / averageHopDistance(i);
            totalHopDistance = totalHopDistance + minHopCount(i, j) / averageHopDistance(i);
        end


        if totalHopDistance > 0

            estimatedUnknownNodePositions_again(j, 1) = sumX / totalHopDistance;
            estimatedUnknownNodePositions_again(j, 2) = sumY / totalHopDistance;
        end
    end
end

    
    
    numAnchors_updated_itr=numAnchors;
    numUnknownNodes_updated_itr=numUnknownNodes;
    anchorPositions_updated_itr=anchorPositions;
    unknownNodePositions_updated_itr=unknownNodePositions;

    updatedAnchorPositions = removeCollinearAnchors(anchorPositions_updated_itr);

anchorPositions_updated_itr=updatedAnchorPositions;
numAnchors_updated_itr=length(updatedAnchorPositions);
function updatedAnchorPositions = removeCollinearAnchors(anchorPositions_updated_itr)
    updatedAnchorPositions = anchorPositions_updated_itr; % Initialize with the original positions
    i = 1;

    while i <= size(updatedAnchorPositions, 1) - 2
        A = updatedAnchorPositions(i, :);
        for j = i + 1:size(updatedAnchorPositions, 1) - 1
            B = updatedAnchorPositions(j, :);
            for k = j + 1:size(updatedAnchorPositions, 1)
                C = updatedAnchorPositions(k, :);
                if isCollinearTriplet(A, B, C)
                    updatedAnchorPositions(k, :) = [];
                    break; 
                end
            end
            if size(updatedAnchorPositions, 1) < size(anchorPositions, 1)
                break; 
            end
        end
        if size(updatedAnchorPositions, 1) < size(anchorPositions, 1)
            continue; 
        end
        i = i + 1;
    end
end

function isCollinear = isCollinearTriplet(A, B, C)
    area = abs(0.5 * (A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2))));
    isCollinear = area < 1e-6;
end




    imp_accuracy = cell(maxIterations, 1);
for iteration = 1:maxIterations %loop run 
   
 RSSI_itr = generateRSSIValues_itr(numAnchors_updated_itr, numUnknownNodes_updated_itr, anchorPositions_updated_itr, unknownNodePositions_updated_itr,communicationRange);
  
  minHopCount_itr = findMinHopCount_itr(RSSI_itr, numAnchors_updated_itr, numUnknownNodes_updated_itr);

   
 
    IMCDS_itr = ImprovedMinimumConnectedDominationSet_itr(numAnchors_updated_itr, numUnknownNodes_updated_itr, anchorPositions_updated_itr, unknownNodePositions_updated_itr, minHopCount_itr, RSSI_itr);
   
    
     averageHopDistance_again_itr = calculateAverageHopDistance_again_itr(IMCDS_itr, unknownNodePositions_updated_itr, minHopCount_itr);


    estimatedUnknownNodePositions_again_itr = calculateUnknownNodeCoordinates_again_itr(IMCDS_itr, minHopCount_itr, averageHopDistance_again_itr);
    
    
    
 
       figure;
    errors_again_itr = sqrt(sum(((unknownNodePositions_updated_itr - estimatedUnknownNodePositions_again_itr)).^2, 2));
    disp('Positional Errors for Each Unknown Node:');
    disp(errors_again_itr);
    hold on;
    plot(1:numUnknownNodes, errors_again, '-.o', 'Color', color1, 'LineWidth', 2, 'MarkerSize', 8);
    plot(1:numUnknownNodes_updated_itr, errors_again_itr, '-*', 'Color', color2, 'LineWidth',2, 'MarkerSize', 8);
  
    title('Accuracy Graph for Each Unknown Node');
    xlabel('Unknown Node Index');
    ylabel('Positioning Error');
    legend('MCDS-DV-HOP', 'AS-IDV-Hop');
    grid on;
  

sizeUpdated = size(unknownNodePositions_updated_itr);
sizeEstimated = size(estimatedUnknownNodePositions_again_itr);

if isequal(sizeUpdated, sizeEstimated)
    differences = abs(unknownNodePositions_updated_itr - estimatedUnknownNodePositions_again_itr);
    
    threshold = 0.9;
    accuratePredictions = differences < threshold;
    numAccuratePredictions = sum(accuratePredictions(:));
    accuracy = numAccuratePredictions / numel(unknownNodePositions_updated_itr);
    accuracy_val= num2str(accuracy * 100);
else
   
    disp('Error: Matrix dimensions must agree.');
end

   
    imp_accuracy{iteration} = accuracy_val;
   
    nodeIndices_itr = (1:size(errors_again_itr, 1))';

    extendedData_itr = [nodeIndices_itr, errors_again_itr];

    sortedData_itr = sortrows(extendedData_itr, 2);       
 
    lowestLocalizationError_itr = sortedData_itr(1:5, 1:2);
    
   lowest_error_list_itr = lowestLocalizationError_itr(:, 1);
   

if exist('unknownNodePositions_updated_itr', 'var')

    validIndices = find(lowest_error_list_itr <= size(unknownNodePositions_updated_itr, 1));

    if ~isempty(validIndices)
  
        removedRows = unknownNodePositions_updated_itr(validIndices, :);
        unknownNodePositions_updated_itr(validIndices, :) = [];
     
        disp(['Removed ', num2str(size(removedRows, 1)), ' rows.']);
    else
        disp('No rows to remove.');
    end
else
    disp('Error: unknownNodePositions_updated_itr.');
   
    break; 
end


   
  
    assistantAnchorNodes_itr = []; 
    assistantAnchorNodes_itr = [assistantAnchorNodes_itr; removedRows];
  

    numAnchors_updated_itr = numAnchors_updated_itr + length(assistantAnchorNodes_itr);

   
   numUnknownNodes_updated_itr=length(unknownNodePositions_updated_itr);
   

   
     disp(numUnknownNodes_updated_itr);
   
     anchorPositions_updated_itr = [anchorPositions_updated_itr; assistantAnchorNodes_itr];
   
   
   
   
   
    
   
   
   
end





numericArray = cellfun(@str2double, imp_accuracy);





end

function RSSI_itr = generateRSSIValues_itr(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions,communicationRange)

    P_tr = 20; 
    P_loss = 30; 
    d_0 = 1; 
    tau = 2; 
    RSSI_itr = zeros(numAnchors, numUnknownNodes);

    for i = 1:numAnchors
        for j = 1:numUnknownNodes
            distance = sqrt((anchorPositions(i, 1) - unknownNodePositions(j, 1))^2 + (anchorPositions(i, 2) - unknownNodePositions(j, 2))^2);
            RSSI_itr(i, j) = P_tr - P_loss - 10 * tau * log10(distance / d_0);
        end
    end
end
 
function minHopCount_itr = findMinHopCount_itr(RSSI_itr, numAnchors_updated_itr, numUnknownNodes_updated_itr)

    minHopCount_itr = zeros(numAnchors_updated_itr, numUnknownNodes_updated_itr);

    for j = 1:numUnknownNodes_updated_itr
        hopCount_itr = inf(numAnchors_updated_itr, 1);

        [~, anchorIdx_itr] = max(RSSI_itr(:, j));
        hopCount_itr(anchorIdx_itr) = 0;

        for i = 1:numAnchors_updated_itr
            for k = 1:numAnchors_updated_itr
                if RSSI_itr(k, j) > RSSI_itr(i, j)
                    hopCount_itr(i) = min(hopCount_itr(i), hopCount_itr(k) + 1);
                end
            end
        end

       
        minHopCount_itr(:, j) = hopCount_itr;
    end
end




function IMCDS_itr = ImprovedMinimumConnectedDominationSet_itr(numAnchors, numUnknownNodes, anchorPositions, unknownNodePositions, minHopCount, RSSI)
    
    IMCDS_Indices = [];

    minRequiredAnchors = ceil(numAnchors / 3);

    connectivity = sum(minHopCount ~= inf, 2); 
    [~, sortedIndices] = sort(connectivity, 'descend'); 

    coveredUnknownNodes = false(numUnknownNodes, 1);
    for i = 1:numAnchors
        anchorIdx = sortedIndices(i);
        connectedUnknowns = minHopCount(anchorIdx, :) ~= inf;

        isNewNodeCovered = any(~coveredUnknownNodes & connectedUnknowns);
        
 
        isBelowMinThreshold = length(IMCDS_Indices) < minRequiredAnchors;

        isNewNodeCovered = isscalar(isNewNodeCovered) && isNewNodeCovered;
        isBelowMinThreshold = isscalar(isBelowMinThreshold) && isBelowMinThreshold;
        if isNewNodeCovered || isBelowMinThreshold
            IMCDS_Indices = [IMCDS_Indices, anchorIdx];
            coveredUnknownNodes = coveredUnknownNodes | connectedUnknowns;
        end
allNodesCovered = isscalar(all(coveredUnknownNodes))&& coveredUnknownNodes;

       
            if((allNodesCovered &&length(IMCDS_Indices))>= minRequiredAnchors)
            break;
        end
    end


    IMCDS_itr = anchorPositions(IMCDS_Indices, :);
end
function averageHopDistance_again_itr = calculateAverageHopDistance_again_itr(anchorPositions, unknownNodePositions, minHopCountRecalculate)
   
    numAnchors = size(anchorPositions, 1);
    disp(minHopCountRecalculate);
    numUnknownNodes = size(unknownNodePositions, 1);

  
    if size(minHopCountRecalculate, 1) < numAnchors || size(minHopCountRecalculate, 2) < numUnknownNodes
        error('minHopCount does not have the correct dimensions.');
    end

    averageHopDistance_again_itr = zeros(numAnchors, 1);

    for i = 1:numAnchors
        totalCorrectedDistance = 0;
        totalCount = 0;

        for j = 1:numUnknownNodes
            if minHopCountRecalculate(i, j) > 0
                actualDistance = sqrt((anchorPositions(i, 1) - unknownNodePositions(j, 1))^2 + (anchorPositions(i, 2) - unknownNodePositions(j, 2))^2);
                correctedDistance = actualDistance / minHopCountRecalculate(i, j); % Correction factor
                totalCorrectedDistance = totalCorrectedDistance + correctedDistance;
                totalCount = totalCount + 1;
            end
        end

        if totalCount > 0
            averageHopDistance_again_itr(i) = totalCorrectedDistance / totalCount;
        end
    end
end

function estimatedUnknownNodePositions_again_itr = calculateUnknownNodeCoordinates_again_itr(anchorPositions, minHopCount, averageHopDistance)
    numUnknownNodes = size(minHopCount, 2);
    estimatedUnknownNodePositions_again_itr = zeros(numUnknownNodes, 2);

    for j = 1:numUnknownNodes
        sumX = 0;
        sumY = 0;
        totalHopDistance = 0;

        for i = 1:size(anchorPositions, 1)
            % Incorporate hop count and average hop distance into the calculation
            sumX = sumX + anchorPositions(i, 1) * minHopCount(i, j) / averageHopDistance(i);
            sumY = sumY + anchorPositions(i, 2) * minHopCount(i, j) / averageHopDistance(i);
            totalHopDistance = totalHopDistance + minHopCount(i, j) / averageHopDistance(i);
        end

        if totalHopDistance > 0
            estimatedUnknownNodePositions_again_itr(j, 1) = sumX / totalHopDistance;
            estimatedUnknownNodePositions_again_itr(j, 2) = sumY / totalHopDistance;
        end
    end
end