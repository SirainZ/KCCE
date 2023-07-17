function result = Ncut(W, label)
    
    l = length(unique(label));
    result = 0;
    for i = 1:l
        Vl = find(label == i);
        r = W(Vl, :);
        V = sum(sum(r));
        r(:, Vl) = 0;
        s = sum(sum(r));
        result = result + s / V;
    end

end