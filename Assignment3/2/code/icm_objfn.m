% Function to compute value of posterior for a given data (noisy) image (given) and
% present estimate of noiseless image (guess)
% The prior_wt is used to weight the prior
% gamma_param is the gamma to be used in case of a Huber or log prior,
% ignored otherwise
% pot_type is the type of potential function used - it can be
%       'quad' for quadratic
%       'huber' for discontinuity-adaptive huber
%        'log' for discontinuity-adaptive  
function val = icm_objfn(given, guess, pot_type, prior_wt, gamma_param)

    if prior_wt<0 || prior_wt>1
        error('Prior weight must be between 0 and 1');
    end
    
    [r,c] = size(given);
    
    lhood = sumsqr(abs(given-guess)); 
    
    up_neigh = abs(guess - circshift(guess,[1,0]));  % difference between pixel and upwards neighbour
    down_neigh = abs(guess - circshift(guess,[r-1,0])); % difference between pixel and downwards neighbour
    left_neigh = abs(guess - circshift(guess,[0,1]));   % difference between pixel and left neighbour
    right_neigh = abs(guess - circshift(guess, [0,c-1]));    % difference between pixel and right neighbour
    
    if strcmp(pot_type,'quad')
        pot_fn = sumsqr(up_neigh)+sumsqr(down_neigh)+sumsqr(left_neigh)+sumsqr(right_neigh); 
        
    elseif strcmp(pot_type, 'huber')
        pot_up = (up_neigh<=gamma_param).*0.5*up_neigh.^2 + (up_neigh>gamma_param).*(gamma_param.*up_neigh -0.5*gamma_param^2);
        pot_down = (down_neigh<=gamma_param).*0.5*down_neigh.^2 + (down_neigh>gamma_param).*(gamma_param.*down_neigh -0.5*gamma_param^2);
        pot_left = (left_neigh<=gamma_param).*0.5*left_neigh.^2 + (left_neigh>gamma_param).*(gamma_param.*left_neigh -0.5*gamma_param^2);
        pot_right = (right_neigh<=gamma_param).*0.5*left_neigh.^2 + (right_neigh>gamma_param).*(gamma_param.*right_neigh -0.5*gamma_param^2);        
        pot_fn = sum(sum(pot_up + pot_down + pot_left + pot_right));
   
    elseif strcmp(pot_type, 'log')
        pot_up = gamma_param*up_neigh - gamma_param^2*log(1+up_neigh/gamma_param); 
        pot_down = gamma_param*down_neigh - gamma_param^2*log(1+down_neigh/gamma_param); 
        pot_left = gamma_param*left_neigh - gamma_param^2*log(1+left_neigh/gamma_param); 
        pot_right = gamma_param*right_neigh - gamma_param^2*log(1+right_neigh/gamma_param); 
        pot_fn = sum(sum(pot_up + pot_down + pot_left + pot_right));

    end
    
    val = (1-prior_wt)*lhood + prior_wt*pot_fn;
    
end