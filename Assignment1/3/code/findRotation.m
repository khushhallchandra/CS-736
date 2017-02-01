function [ R ] = findRotation(img1, img2)
% This function returns the rotation matrix
    [U,~,V] = svd(img1*img2');
    R = V*U';
    if(det(R) == -1)
        R = V*[1 0;0 -1]*U';
    end
end