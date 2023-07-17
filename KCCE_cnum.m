function v = KCCE_cnum(X, sigma, kmin, interval, kmax, Pnum)

[N, ~] = size(X);
[ind_kmeans, C] = kmeans(X, Pnum);

[num, d] = size(C);
D = squareform(pdist(C))  ;
W = exp(-(D / sigma).^2);
%     W = 1 - D;
T = 1./sum(W);
W = sqrt(T)' .* W .* sqrt(T);


klist = kmin: interval: kmax;
cnum = zeros([length(klist), 1]);
cnum(1) = num;

Wk = W ^ klist(1);
Ws = W ^ interval;
for i = 2:length(klist)
    Wk = Wk * Ws;
    [~, idx] = max(Wk, [], 1);
    c = find((1:num) == idx)';
    cnum(i) = length(c);
end

v = [klist', cnum];

end


