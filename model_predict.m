function yhat = model_predict(u1, u2, params)


M = size(params,1);
N = numel(u1);
yhat = zeros(N,1);

for i = 1:M
    w  = params(i,1);
    c1 = params(i,2);
    c2 = params(i,3);
    s1 = params(i,4);
    s2 = params(i,5);

    G = exp(-(((u1-c1).^2)./(2*s1^2) + ((u2-c2).^2)./(2*s2^2)));

    yhat = yhat + w*G;
end
end
