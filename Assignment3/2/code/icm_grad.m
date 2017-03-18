% Function to compute gradient of posterior for given data (noisy) image
% (given) and current estimate of noiseless image (guess)
% The prior_wt is used to weight the prior
% gamma_param is the gamma to be used in case of a Huber or log prior,
% ignored otherwise
% pot_type is the type of potential function used - it can be
%       'quad' for quadratic
%       'huber' for discontinuity-adaptive huber
%        'log' for discontinuity-adaptive  

function grad = icm_grad(given, guess, pot_type, prior_wt, gamma_param)

    if prior_wt<0 || prior_wt>1
        error('Prior weight must be between 0 and 1');
    end
    
    % Vectorizing images
    given_vec = given(:);
    guess_vec = guess(:);
    
    [r,c] = size(given);
    
    deriv_lhood = 2*(guess - given);    % Derivative of likelihood term

    up_neigh = guess - circshift(guess,[1,0]);  % difference between pixel and upwards neighbour
    down_neigh = guess - circshift(guess,[r-1,0]); % difference between pixel and downwards neighbour
    left_neigh = guess - circshift(guess,[0,1]);   % difference between pixel and left neighbour
    right_neigh = guess - circshift(guess, [0,c-1]);    % difference between pixel and right neighbour
    
    % deriv_pot is the derivative of the potential function
    if strcmp(pot_type, 'quad')
        deriv_pot = 2*(up_neigh + down_neigh + left_neigh + right_neigh);    
    
    elseif strcmp(pot_type, 'huber')
        deriv_up = (abs(up_neigh)<=gamma_param).*up_neigh + (abs(up_neigh)>gamma_param).*gamma_param.*sign(up_neigh);
        deriv_down = (abs(down_neigh)<=gamma_param).*down_neigh + (abs(down_neigh)>gamma_param).*gamma_param.*sign(down_neigh);
        deriv_left = (abs(left_neigh)<=gamma_param).*left_neigh + (abs(left_neigh)>gamma_param).*gamma_param.*sign(left_neigh);
        deriv_right = (abs(right_neigh)<=gamma_param).*right_neigh + (abs(right_neigh)>gamma_param).*gamma_param.*sign(right_neigh);
        deriv_pot = deriv_up + deriv_down + deriv_left + deriv_right;
    
    elseif strcmp(pot_type, 'log')
        deriv_up = gamma_param.*up_neigh./(gamma_param+abs(up_neigh));
        deriv_down = gamma_param.*down_neigh./(gamma_param+abs(down_neigh));
        deriv_left = gamma_param.*left_neigh./(gamma_param+abs(left_neigh));
        deriv_right = gamma_param.*right_neigh./(gamma_param+abs(right_neigh));
        deriv_pot = deriv_up + deriv_down + deriv_left + deriv_right;
        
    end
        
    grad = (1-prior_wt)*deriv_lhood + prior_wt*deriv_pot;

end