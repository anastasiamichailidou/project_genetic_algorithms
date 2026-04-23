function idx = tournament_select(fitness, k)


P = numel(fitness);

candidates = randi(P, [k,1]);   

[~, bestPos] = min(fitness(candidates));
idx = candidates(bestPos);       
end
