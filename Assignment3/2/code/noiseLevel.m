function sigma = getNoiseLevel(imageNoisy)
p1 = imageNoisy(1:130, 1:40);
p2 = imageNoisy(169:256, 1:30);
p3 = imageNoisy(200:256, 230:256);
p4 = imageNoisy(1:30, 230:256);
% Make a row vector
patch = [reshape(p1, 1, []) reshape(p2, 1, []) reshape(p3, 1, []) reshape(p4, 1, []);];

sigma = std(real(patch));
end
