function [u] = F(m,s)

u = (1+s.^(-m)).^-1;
if m == 0
    u =1;
end
end

