% function to compute Radon transform

function radon_trans = myRadonTrans(img, t, theta, step_size)
    
    radon_trans = zeros(length(t),length(theta));
    
    for i=1:length(t)
        for j=1:length(theta)
            radon_trans(i,j) = myIntegration(img, t(i), theta(j), step_size);
        end
    end

end