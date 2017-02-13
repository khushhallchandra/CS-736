% function to compute line integral along a line specified by normal
% distance t and with the normal at an angle theta degrees

function integral = myIntegration(img, t, theta, step_size)

    [ver_size, hor_size] = size(img);

    centre_x = floor(hor_size/2);
    centre_y = floor(ver_size/2);
    
    if (t*cosd(theta)<(1-centre_x)||t*cosd(theta)>centre_x||t*sind(theta)<(1-centre_x)||t*sind(theta)>centre_y) % Line doesn't intersect image
        integral = 0;
    else
        % Determining limits of integration
        s_poss = floor(sqrt(centre_x^2+centre_y^2)/2);
        s_all = [-1*s_poss:step_size:s_poss];
        x_all = t*cosd(theta) - s_all*sind(theta);
        y_all = t*sind(theta) + s_all*cosd(theta);

        mask = (x_all>=(1-centre_x)).*(x_all<=centre_x).*(y_all>=(1-centre_y)).*(y_all<=centre_y);     % Ones indicate which values of s, x and y lie within the image

        %s = nonzeros(s_all.*mask);
        ind_min = find(mask,1,'first');
        ind_max = find(mask,1,'last');
        x_s = x_all(ind_min:ind_max);
        y_s = y_all(ind_min:ind_max);
        %fprintf('t = %d, theta = %d, length_x = %d, length_y = %d\n', t, theta, length(x_s), length(y_s));

        % Linear interpolation
        [X,Y] = meshgrid(-centre_x+1:centre_x, -centre_y+1:centre_y);
        img_inter = interp2(X,Y,img,x_s,y_s,'linear');
        integral = sum(img_inter)*step_size;
    end
%     % Determing the limits of s for the line integral
%     integral_zero = false;    % Zero if line does not cut image 
%     
%     if(theta==0)
%         s_max = centre_y;
%         s_min = 1-centre_y;
%         if(t>centre_x || t<1-centre_x) 
%             integral_zero = true;
%         end
% 
%     elseif(theta>0 && theta<90)
%         if(t>=0)
%             s_min = (t*cosd(theta)-centre_x)/sind(theta);
%             s_max = (centre_y-t*sind(theta))/cosd(theta);
%         else
%             s_max = (t*cosd(theta)-(-centre_x+1))/sind(theta);
%             s_min = (-centre_y+1-t*sind(theta))/cosd(theta);
%         end
%     elseif(theta==90)
%         s_max = centre_x;
%         s_min = 1-centre_x;
%         if(t>centre_y || t<1-centre_y) 
%             integral_zero = true;
%         end
%     else
%         if(t>=0)
%             s_min = (centre_y-t*sind(theta))/cosd(theta);
%             s_max = (t*cosd(theta)-(-centre_x+1))/sind(theta);
%         else
%             s_max = ((-centre_y+1)-t*sind(theta))/cosd(theta);
%             s_min = (t*cosd(theta)-centre_x)/sind(theta);
%         end
%     end
%     
%     % Indicates that line is not cutting image
%     if(s_max<s_min) 
%         integral_zero = true;
%     end
% 
%     % Linear interpolation
%     s = [s_min:step_size:s_max];
%     x_s = t*cosd(theta)-s*sind(theta);
%     y_s = t*sind(theta)+s*cosd(theta);
% 
%     [X_S, Y_S] = meshgrid(x_s, y_s);
% 
%     [X,Y] = meshgrid(-centre_x+1:centre_x, -centre_y+1:centre_y);
%     img_inter = interp2(X,Y,img,x_s,y_s,'linear');
% 
%     integral = sum(img_inter)*step_size;
%     if(integral_zero)
%         integral = 0;
end

