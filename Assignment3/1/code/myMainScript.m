%% Q1
close all
clear all
clc

%% Loading data
img_data = load('../data/assignmentImageDenoisingPhantom.mat');
img_orig = img_data.imageNoiseless;
img_noisy = img_data.imageNoisy;

%% (a) RRMSE between noisy and noiseless images
rrmse_noise = rrmse_calc(img_noisy, img_orig);
fprintf('The RRMSE between the noisy and noiseless images is %.4f\n', rrmse_noise);

%% Gradient descent

%% Quadratic prior
%% Optimal alpha
alpha_opt = 0.21;            % Weight given to prior
alpha = alpha_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'quad',alpha,1); 
obj_quad = [val_curr];          % Array storing value 

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
    
    obj_quad = [obj_quad val_curr];
end
img_denoised_quad = img_curr;
rrmse_denoised = rrmse_calc(img_denoised_quad, img_orig);

%% Nearby parameter values

alpha = 0.8*alpha_opt;
step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'quad',alpha,1); 

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
end
img_denoised1 = img_curr;
rrmse_denoised1 = rrmse_calc(img_denoised1, img_orig);


alpha = 1.2*alpha_opt;
step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'quad',alpha,1); 

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
end
img_denoised2 = img_curr;
rrmse_denoised2 = rrmse_calc(img_denoised2, img_orig);

fprintf('The optimal alpha is %.3f with RRMSE = %.5f\n', alpha_opt, rrmse_denoised);
fprintf('The RRMSE for alpha = %.3f is %.5f\n', 0.8*alpha_opt, rrmse_denoised1);
fprintf('The RRMSE for alpha = %.3f is %.5f\n', 1.2*alpha_opt, rrmse_denoised2);


%% Huber prior
%% Optimal alpha and gamma
alpha_opt = 0.995;            % Weight given to prior
gamma_opt = 0.0052;
alpha = alpha_opt;
gamma_param = gamma_opt;


step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
obj_huber = [val_curr];          % Array storing value 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
    obj_huber = [obj_huber val_curr];
end
img_denoised_huber = img_curr;
rrmse_denoised = rrmse_calc(img_denoised_huber, img_orig);
fprintf('RRMSE for optimal alpha = %.4f, optimal gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised)

%% Nearby parameter values
alpha = 0.8*alpha_opt;
gamma_param = gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised10 = img_curr;
rrmse_denoised10 = rrmse_calc(img_denoised10, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised10)


alpha = min(1.2*alpha_opt,1);
gamma_param = gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised20 = img_curr;
rrmse_denoised20 = rrmse_calc(img_denoised20, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised20)


alpha = alpha_opt;
gamma_param = 0.8*gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised01 = img_curr;
rrmse_denoised01 = rrmse_calc(img_denoised01, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised01)


alpha = alpha_opt;
gamma_param = 1.2*gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'huber',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'huber',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'huber',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised02 = img_curr;
rrmse_denoised02 = rrmse_calc(img_denoised02, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised02)


%% Log prior
%% Optimal alpha and gamma
alpha_opt = 0.978;            % Weight given to prior
gamma_opt = 0.005;
alpha = alpha_opt;
gamma_param = gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
obj_log = [val_curr];          % Array storing value 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
    obj_log = [obj_log val_curr];
end
img_denoised_log = img_curr;
rrmse_denoised = rrmse_calc(img_denoised_log, img_orig);
fprintf('RRMSE for optimal alpha = %.4f, optimal gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised)



%% Nearby parameter values
alpha = 0.8*alpha_opt;
gamma_param = gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised10 = img_curr;
rrmse_denoised10 = rrmse_calc(img_denoised10, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised10);


alpha = min(1.2*alpha_opt,1);
gamma_param = gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised20 = img_curr;
rrmse_denoised20 = rrmse_calc(img_denoised20, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised20);


alpha = alpha_opt;
gamma_param = 0.8*gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised01 = img_curr;
rrmse_denoised01 = rrmse_calc(img_denoised01, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised01);



alpha = alpha_opt;
gamma_param = 1.2*gamma_opt;

step = 1;  % initial step size
img_curr = img_noisy;   % present guess for denoised image
stop = 1e-6;  % Stopping criteria on step size    
val_curr = icm_objfn(img_noisy,img_curr,'log',alpha,gamma_param); 
iter = 0;
while(iter<=20)    
    img_grad = icm_grad(img_noisy,img_curr,'log',alpha, gamma_param);
    img_upd = img_curr - step*img_grad;     % Image to be tested
    val_upd = icm_objfn(img_noisy,img_upd,'log',alpha,gamma_param);
%     fprintf('alpha = %.3f, gamma = %.3f, Value of current obj fn = %.5f, possible obj fn = %.5f\n', alpha, gamma_param, val_curr, val_upd);
    
    if(val_upd<val_curr)            % if objective function decreases in values
        img_curr = img_upd;         % update estimate of image
        val_curr = val_upd;         % Update value of objective function
        step = 1.1*step;            % Increase step size by 10%
        
    else
        step = 0.5*step;            % Decrease step size by 50% if objective function does not decrease
    end
    iter = iter+1;
end
img_denoised02 = img_curr;
rrmse_denoised02 = rrmse_calc(img_denoised02, img_orig);
fprintf('RRMSE for alpha = %.4f, gamma = %.4f is %.5f\n', alpha, gamma_param, rrmse_denoised02);

%% Displaying results
range_max = max([max(abs(img_orig(:))), max(abs(img_noisy(:))), max(abs(img_denoised_quad(:))), max(abs(img_denoised_huber(:))), max(abs(img_denoised_log(:)))]);
range_min = min([min(abs(img_orig(:))), min(abs(img_noisy(:))), min(abs(img_denoised_quad(:))), min(abs(img_denoised_huber(:))), min(abs(img_denoised_log(:)))]);

figure(1);
imshow(abs(img_orig),[range_min, range_max]);
title('Noiseless image');
colorbar

figure(2);
imshow(abs(img_noisy),[range_min, range_max]);
title('Noisy image');
colorbar

figure(3);
imshow(abs(img_denoised_quad),[range_min, range_max]);
title('Denoised image using quadratic prior');
colorbar

figure(4);
imshow(abs(img_denoised_huber),[range_min, range_max]);
title('Denoised image using Huber prior');
colorbar

figure(5);
imshow(abs(img_denoised_log),[range_min, range_max]);
title('Denoised image using log prior');
colorbar

%% Plots of objective function values
%% For quadratic prior
figure(6);
plot(obj_quad)
xlabel('Number of iterations')
ylabel('Value of objective function');
title('For quadratic prior');

%% For Huber prior
figure(7);
plot(obj_huber)
xlabel('Number of iterations')
ylabel('Value of objective function');
title('For Huber prior');

%% For log prior
figure(8);
plot(obj_log)
xlabel('Number of iterations')
ylabel('Value of objective function');
title('For log prior');