clear all
close all
clc
%% Load data

filename = '../data/ellipses2D.mat';
data = load(filename);
numPoints = data.numOfPoints;
numSets = data.numOfPointSets;
pointSets = data.pointSets;

%% Plot initial point set
figure(1);
hold on;
title('Initial Point Sets');
for i=1:numSets
    scatter(pointSets(1,:,i),pointSets(2,:,i),10,rand(1,3));
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
title('Plot of variances');
plot(eigvals);
hold off;

%% Principal modes of shape variation
sd2 = sqrt(eigvals(2));
sd1 = sqrt(eigvals(1));

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
hold on;
plot(pmv11(1,:),pmv11(2,:),'--b*');
plot(pmv12(1,:),pmv12(2,:),'--b*');
plot(meanShape(1,:),meanShape(2,:),'--r*');
title('1st mode of variation');
hold off;

%% Plot 2nd mode of variation

figure(5)
hold on;
plot(pmv21(1,:),pmv21(2,:),'--b*');
plot(pmv22(1,:),pmv22(2,:),'--b*');
plot(meanShape(1,:),meanShape(2,:),'--r*');
title('2nd mode of variation');
hold off;
