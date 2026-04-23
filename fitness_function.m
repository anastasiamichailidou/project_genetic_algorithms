function J = fitness_function(chrom, M, u1, u2, y)


params = decode_chrom(chrom, M);
yhat = model_predict(u1, u2, params);

J = mean((y - yhat).^2);
end
