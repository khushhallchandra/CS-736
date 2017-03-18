%% Q2
close all
clear all
clc

%% Loading data
img_data = load('../data/assignmentImageDenoisingBrainNoisy.mat');
img_noisy = img_data.imageNoisy;
% a) Estimate the noise level
% noise level
sd = noiseLevel(img_noisy);
display(sd);


%% Gradient descent

%% Quadratic prior
alpha_opt = 0.5;            % Weight given to prior
alpha = alpha_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'quad',alpha,1); 
obj_quad1 = [val_curr];          % Array storing value 

while(step>stop)    
    img_grad = icm_grad(img_noisy,img_curr,'quad',alpha,1);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'quad',alpha,1);
%     fprintf('alpha = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    
    obj_quad1 = [obj_quad1 val_curr];
end
img_denoised1 = img_curr;
% Optimum alpha
display(alpha_opt);

%% Huber prior
alpha_opt = 0.4;            % Weight given to prior
alpha = alpha_opt;
gamma_opt = 0.3;
gamma_param = gamma_opt;
%rrmse_min = rrmse_noise;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image   
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
obj_quad2 = [val_curr];          % Array storing value 
iter = 0;

while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter + 1;
    obj_quad2 = [obj_quad2 val_curr];
end
img_denoised2 = img_curr;

% Optimum alpha
display(alpha_opt);
% Optimum gamma
display(gamma_opt);

%% Discontinuity-adaptive prior
alpha_opt = 0.4;            % Weight given to prior
alpha = alpha_opt;
gamma_opt = 0.5;
gamma_param = gamma_opt;
%rrmse_min = rrmse_noise;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-5;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
obj_quad3 = [val_curr];          % Array storing value 
iter = 0;

while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter + 1;
    obj_quad3 = [obj_quad3 val_curr];
end
img_denoised3 = img_curr;
% Optimum alpha
display(alpha_opt);
% Optimum gamma
display(gamma_opt);


figure(1)
h(1) = subplot(1,4,1);      
imshow(abs(img_noisy),[0,1]);
title('Noisy image');

h(2) = subplot(1,4,2); 
imshow(abs(img_denoised1),[0,1]);
title('Quadratic prior');

h(3) = subplot(1,4,3);
imshow(abs(img_denoised2),[0,1]);
title('Huber Prior');

h(4) = subplot(1,4,4);
imshow(abs(img_denoised3),[0,1]);
title('discontinuity-adaptive');

suptitle('Noisy and Denoised Images');
pause(5);

display('All the four images are on same colormap');

figure(2)

h(1) = subplot(1,3,1); 
plot(obj_quad1);
title('Quadratic prior');

h(2) = subplot(1,3,2);
plot(obj_quad2);
title('Huber Prior');

h(3) = subplot(1,3,3);
plot(obj_quad3);
title('discontinuity-adaptive');

suptitle('Objective-function values');
