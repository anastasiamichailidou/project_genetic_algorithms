function chrom = mutate_gaussian(chrom, pmut, sigma_mut)


mask = rand(size(chrom)) < pmut;
chrom(mask) = chrom(mask) + sigma_mut*randn(1, sum(mask));

end
