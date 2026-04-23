function child = crossover_blend(p1, p2, alpha)

beta = (1+2*alpha)*rand(size(p1)) - alpha;  
child = (1-beta).*p1 + beta.*p2;
end
