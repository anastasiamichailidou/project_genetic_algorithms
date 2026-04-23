function chrom = repair_bounds(chrom, M, u1_min, u1_max, u2_min, u2_max)
params = decode_chrom(chrom, M);


params(:,2) = min(max(params(:,2), u1_min), u1_max);
params(:,3) = min(max(params(:,3), u2_min), u2_max);

sig_min = 0.05;
sig_max = 2.0;
params(:,4) = min(max(abs(params(:,4)), sig_min), sig_max);
params(:,5) = min(max(abs(params(:,5)), sig_min), sig_max);

chrom = reshape(params', [1, 5*M]);
end
