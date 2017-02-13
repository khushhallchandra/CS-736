function p = myFilter( radon_trans,filter, L)

n = size(radon_trans,1);
fft_len = 1024;

filt = 2*( 0:(fft_len/2) )./fft_len;
w = 2*pi*(0:size(filt,2)-1)/fft_len;   

if(strcmp(filter, 'Ram-Lak'))
    %
elseif(strcmp(filter, 'Shepp-Logan'))
   filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*L))./(w(2:end)/(2*L)));

elseif(strcmp(filter, 'Cosine'))
   filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*L));

end

% Windowing the filter
filt(w>pi*L) = 0;                      

% Extending filter to -ve frequencies
filt = [filt' ; filt(end-1:-1:2)'];    

% Filtering
p = fft(radon_trans,fft_len); 

for i = 1:size(p,2)
   p(:,i) = p(:,i).*filt; 
end

p = real(ifft(p));     % p is the filtered projections
p(n+1:end,:) = [];   % Truncate the filtered projections
