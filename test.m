clc
clear
N = 10000;

propensities = [ ...
    183   4   0   4   2   2   0   4   4   9   2   1   2   0  15   3   2   3   3 ...
    15   3   1   5   3   9   3   6   9   9 152  10   7   6   5  10   7   2   4 ...
    11  49   6   6   5  12 108   7   6   4   7  49   3  10   1   5  23   3   8 ...
    15   7 348 ...
];



sumProp = sum(propensities);
probs = propensities / sumProp;

L = length(propensities);
result = zeros(L,1);

all = zeros(L, 1);
for i=1:N
    r = randsample(1:L, 1, true, probs);
    result(r) = result(r) + 1;
    all(i) = r;
end

mean(all)
std(all)

result ./ sum(result)

mean(result)

return

result = zeros(4,1);
for i=1:N
    r = randsample(1:4, 1, true, [0.1, 0.2, 0.3, 0.4]);
    result(r) = result(r) + 1;
end

result ./ sum(result)

mean(result)