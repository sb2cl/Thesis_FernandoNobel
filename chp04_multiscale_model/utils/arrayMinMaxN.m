function [out] = arrayMinMaxN(limits,n)
%% This functions generates a n-array of values between valueMin and valueMax.

valueMin = limits(1);
valueMax = limits(2);

if n == 1
    out = valueMin*(1) + valueMax*(0);
else
    out = valueMin + (0:n-1)*(valueMax - valueMin)/(n-1);
end

end

