function [out] = mergeStruct(A,B)
out = B;

f = fieldnames(A);

for i = 1:length(f)
    out.(f{i}) = A.(f{i});
end

end

