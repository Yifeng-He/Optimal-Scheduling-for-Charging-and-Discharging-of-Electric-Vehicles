function x = shuffle(x)
%SHUFFLE  Shuffle an array of values.
%   x = shuffle(x) shuffles the vector x in random order.

%   Since this function changes the vector x,
%   it creates a copy of it.

n = length(x);
disp n
for index = 1:n

    % a random number between i and n
    r = index + floor((n-index+1) * rand());
    
    % swap elements index and r
    x([index r]) = x([r index]);
end
