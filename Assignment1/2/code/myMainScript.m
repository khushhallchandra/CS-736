clear all
close all
clc

%% Load data

filename = '../data/hands2D.mat';
data = load(filename);
pointSets = data.shapes;
[x, numPoints, numSets] = size(pointSets);

%% Plot initial point set
figure(1);
hold on;
title('Initial Point Sets');
for i=1:numSets
    scatter(pointSets(1,:,i),pointSets(2,:,i),12,rand(1,3));
end
hold off;

%% Preshape space
centroid = sum(pointSets,2)/numPoints;
% subtract centroid 
preshapePointSets = pointSets - repmat(centroid,[1,numPoints,1]);
% scaling 
norms = sqrt(sum(sum(preshapePointSets.^2,2),1)); 
preshapePointSets = preshapePointSets./repmat(norms,[2,numPoints,1]);

%% Compute mean shape

% parametrs 
threshold = 1e-7;

% Initial values
% Let's take first pointSet for initialising meanShape 
meanShape = preshapePointSets(:,:,1);

error = 1;
while(threshold < error)
    for i = 1:numSets
        % Given meanShape, find optimal transformations
        R = findRotation(preshapePointSets(:,:,i), meanShape);
        % Apply the transformation
        preshapePointSets(:,:,i) = R * preshapePointSets(:,:,i);        
    end
    
    % Average all (aligned) pointsets
    optMeanShape = sum(preshapePointSets,3)/numSets;
    % normalize the optMeanShape
    optMeanShape = optMeanShape./sqrt(sum(sum(optMeanShape.^2)));
    error = sqrt(sum(sum((optMeanShape - meanShape).^2)));
    % update prevMeanShape
    meanShape = optMeanShape;
end


%% Plot of computed shape mean, together with all the aligned pointsets
figure(2);
hold on;
title('Mean shape amd aligned point sets');
for i=1:numSets
    scatter(preshapePointSets(1,:,i),preshapePointSets(2,:,i),12,rand(1,3));
end
plot(meanShape(1,:),meanShape(2,:),'--ro','LineWidth',2);
hold off;

%% Covariance
pointSetsNew = bsxfun(@minus, preshapePointSets, meanShape);
vectorizedPointSets = reshape(pointSetsNew, 2*numPoints, numSets);
covariance = vectorizedPointSets*vectorizedPointSets'/numSets;
[V,D] = eig(covariance);
% To get descending eigenvalues
V = fliplr(V);
eigvals = flipud(diag(D));

%% Plot variance
figure(3);
hold on;
title('Plot of variance');
plot(eigvals);
hold off;

%% Principal modes of shape variation
sd3 = sqrt(eigvals(3));
sd2 = sqrt(eigvals(2));
sd1 = sqrt(eigvals(1));

% 2*sd along 3rd eigenvector
pmv31 = meanShape + 2*sd3*reshape(V(:,3),2,numPoints);
% -2*sd along 3rd eigenvector
pmv32 = meanShape - 2*sd3*reshape(V(:,3),2,numPoints);

% 2*sd along 2nd eigenvector
pmv21 = meanShape + 2*sd2*reshape(V(:,2),2,numPoints);
% -2*sd along 2nd eigenvector
pmv22 = meanShape - 2*sd2*reshape(V(:,2),2,numPoints);

% 2*sd along ist eigenvector
pmv11 = meanShape + 2*sd1*reshape(V(:,1),2,numPoints);
% -2*sd along 1st eigenvector
pmv12 = meanShape - 2*sd1*reshape(V(:,1),2,numPoints);

%% Plot 1st mode of variation

figure(4)
h(1) = subplot(1,3,1);      
plot(pmv11(2,:), pmv11(1,:),'--b*');
title('+2sd');

h(2) = subplot(1,3,2);
title('mean shape'); 
plot(meanShape(2,:),meanShape(1,:),'--r*');

h(3) = subplot(1,3,3);
title('-2sd');
plot(pmv12(2,:),pmv12(1,:),'--b*');

linkaxes(h)
xlim([-0.23 0.2])
suptitle('1st mode of variation (90 degree rotated)');
pause(5)

%% Plot 2nd mode of variation

figure(5)
h(1) = subplot(1,3,1);      
plot(pmv21(2,:), pmv21(1,:),'--b*');
title('+2sd');

h(2) = subplot(1,3,2);
title('mean shape'); 
plot(meanShape(2,:),meanShape(1,:),'--r*');

h(3) = subplot(1,3,3);
title('-2sd');
plot(pmv22(2,:),pmv22(1,:),'--b*');

linkaxes(h)
xlim([-0.23 0.2])
suptitle('2nd mode of variation (90 degree rotated)');
pause(5)

%% Plot 3rd mode of variation

figure(6)
h(1) = subplot(1,3,1);      
plot(pmv31(2,:), pmv31(1,:),'--b*');
title('+2sd');

h(2) = subplot(1,3,2);
title('mean shape'); 
plot(meanShape(2,:),meanShape(1,:),'--r*');

h(3) = subplot(1,3,3);
title('-2sd');
plot(pmv32(2,:),pmv32(1,:),'--b*');

linkaxes(h)
xlim([-0.23 0.2])
suptitle('3rd mode of variation (90 degree rotated)');
pause(5)