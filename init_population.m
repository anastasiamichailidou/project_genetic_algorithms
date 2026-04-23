function pop = init_population(P, M, u1_min, u1_max, u2_min, u2_max)


D = 5*M;
pop = zeros(P, D);

for p = 1:P
    params = zeros(M,5);
    params(:,1) = randn(M,1); 
    params(:,2) = u1_min + (u1_max-u1_min)*rand(M,1);
    params(:,3) = u2_min + (u2_max-u2_min)*rand(M,1); 
    params(:,4) = 0.2 + 0.8*rand(M,1);
    params(:,5) = 0.2 + 0.8*rand(M,1);

    pop(p,:) = reshape(params', [1, D]);
end
end
