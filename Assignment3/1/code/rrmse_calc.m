% Function to calculate the RRMSE between ref and noisy 
function rrmse = rrmse_calc(noisy, ref)
    diff = abs(ref)-abs(noisy);
    rrmse = sqrt(sumsqr(diff)/sumsqr(abs(ref)));
end