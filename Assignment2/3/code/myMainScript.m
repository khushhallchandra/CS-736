%% Q3

clear all
close all

%% General comments:
% We would take that contiguous set of 150 angles which captures most of the
% non-zero part of the image. The exact set is difficult to determine from
% observation but can be found through calculation of RRMSE for different
% such sets.

%% CT_Chest data
struct_load = load('../data/CT_Chest.mat');
chest = struct_load.imageAC;

figure(1);
imshow(chest,[]);
title('CT Chest image');

% radonTransform = radon(chest,1:180);
% considering only 150 thetas
thetas = 1:150;
RRMSE = zeros(180,1);

normVal = sqrt(sumsqr(chest));
%Check for all thetas
for i = 1:180

    radonTransform = radon(chest,mod(thetas+i,180));
    inverseRadon = iradon(radonTransform, mod(thetas+i,180), 'linear', 'Ram-Lak', 1,size(chest,1));
    RRMSE(i) = sqrt( sumsqr(inverseRadon - chest) )/normVal;
end

[~,minIndex] = min(RRMSE);
minRadonTransform = radon(chest,thetas+minIndex);
minInverseRadon = iradon(minRadonTransform, mod(thetas+minIndex,180), 'linear', 'Ram-Lak', 1,size(chest,1));
figure(2)
imshow(minInverseRadon,[])
title('Reconstructed image for least RRMSE')
fprintf(' value of theta for least RRMSE is %d \n ',minIndex)
pause(10)
figure(3)
plot(RRMSE)
title('RRMSE for different \thetas ')
xlabel('\theta')
ylabel('RRMSE')
pause(10)
%% myPhantom data
struct_load = load('../data/myPhantom.mat');
myPhantom = struct_load.imageAC;

figure(4);
imshow(myPhantom,[]);
colorbar
title('Phantom image');

%radonTransform = radon(myPhantom,1:180);
% considering only 150 thetas
thetas = 1:150;
RRMSE = zeros(180,1);

normVal = sqrt(sumsqr(myPhantom));
%Check for all thetas
for i = 1:180

    radonTransform = radon(myPhantom,mod(thetas+i,180));
    inverseRadon = iradon(radonTransform, mod(thetas+i,180), 'linear', 'Ram-Lak', 1,size(myPhantom,1));
    RRMSE(i) = sqrt( sumsqr(inverseRadon - myPhantom) )/normVal;
end

[~,minIndex] = min(RRMSE);
minRadonTransform = radon(myPhantom,thetas+minIndex);
minInverseRadon = iradon(minRadonTransform, mod(thetas+minIndex,180), 'linear', 'Ram-Lak', 1,size(myPhantom,1));
figure(5)
imshow(minInverseRadon,[])
colorbar
title('Reconstructed image for least RRMSE')
fprintf(' value of theta for least RRMSE is %d \n ',minIndex)
pause(10)
figure(6)
plot(RRMSE)
title('RRMSE for different \thetas ')
xlabel('\theta')
ylabel('RRMSE')
pause(10)