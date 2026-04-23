clear; clc; close all;

rng(42);  

M_values = [5, 8, 10, 12, 15]; 
num_runs = length(M_values);

P = 50;         
gens = 100;     
k = 3;         
alpha = 0.5;    
pmut = 0.10;    
sigma_mut = 0.10;
elite = 2;

Ntrain = 800;
Ntest  = 400;

u1_min = -1; u1_max = 2;
u2_min = -2; u2_max = 1;

u1_tr = u1_min + (u1_max-u1_min)*rand(Ntrain,1);
u2_tr = u2_min + (u2_max-u2_min)*rand(Ntrain,1);
u1_te = u1_min + (u1_max-u1_min)*rand(Ntest,1);
u2_te = u2_min + (u2_max-u2_min)*rand(Ntest,1);

f = @(u1,u2) sin(u1 + u2).*sin(u2.^2);
y_tr = f(u1_tr, u2_tr);
y_te = f(u1_te, u2_te);

results = struct();
results.M = M_values;
results.train_mse = zeros(num_runs, 1);
results.test_mse = zeros(num_runs, 1);
results.best_params = cell(num_runs, 1);
results.best_hist = cell(num_runs, 1);

for idx = 1:num_runs
    M = M_values(idx);
    fprintf('\n========================================\n');
    fprintf('RUNNING: M = %d (experiment %d/%d)\n', M, idx, num_runs);
    fprintf('========================================\n');
    

    pop = init_population(P, M, u1_min, u1_max, u2_min, u2_max);
    best_hist = zeros(gens,1);
    mean_hist = zeros(gens,1);

    for g = 1:gens
       
        fitness = zeros(P,1);
        for p = 1:P
            fitness(p) = fitness_function(pop(p,:), M, u1_tr, u2_tr, y_tr);
        end
        
      
        [fitness_sorted, sort_idx] = sort(fitness);
        pop_sorted = pop(sort_idx,:);
        
        best_hist(g) = fitness_sorted(1);
        mean_hist(g) = mean(fitness_sorted);
        
  
        newpop = zeros(size(pop));
        newpop(1:elite,:) = pop_sorted(1:elite,:);
        
      
        for p = (elite+1):P
            idx1 = tournament_select(fitness_sorted, k);
            idx2 = tournament_select(fitness_sorted, k);
            parent1 = pop_sorted(idx1, :);
            parent2 = pop_sorted(idx2, :);
            
            child = crossover_blend(parent1, parent2, alpha);
            child = mutate_gaussian(child, pmut, sigma_mut);
            child = repair_bounds(child, M, u1_min, u1_max, u2_min, u2_max);
            newpop(p,:) = child;
        end
        
        pop = newpop;
        
        if mod(g,20)==0
            fprintf("  Gen %d: best=%.6f  mean=%.6f\n", g, best_hist(g), mean_hist(g));
        end
    end
    
 
    fitness = zeros(P,1);
    for p = 1:P
        fitness(p) = fitness_function(pop(p,:), M, u1_tr, u2_tr, y_tr);
    end
    [~, bestIdx] = min(fitness);
    bestChrom = pop(bestIdx,:);
    bestParams = decode_chrom(bestChrom, M);
    

    yhat_tr = model_predict(u1_tr, u2_tr, bestParams);
    yhat_te = model_predict(u1_te, u2_te, bestParams);
    mse_tr = mean((y_tr - yhat_tr).^2);
    mse_te = mean((y_te - yhat_te).^2);
  
    results.train_mse(idx) = mse_tr;
    results.test_mse(idx) = mse_te;
    results.best_params{idx} = bestParams;
    results.best_hist{idx} = best_hist;
    
    fprintf('\n  RESULTS: Train MSE=%.6f | Test MSE=%.6f\n', mse_tr, mse_te);
end

fprintf('\n\n');
fprintf('========================================\n');
fprintf('         ΣΥΓΚΕΝΤΡΩΤΙΚΑ ΑΠΟΤΕΛΕΣΜΑΤΑ     \n');
fprintf('========================================\n');
fprintf('  M  | Train MSE  | Test MSE   | Ratio\n');
fprintf('-----|------------|------------|---------\n');
for idx = 1:num_runs
    ratio = results.test_mse(idx) / results.train_mse(idx);
    fprintf(' %2d  | %.6f | %.6f | %.3f\n', ...
        results.M(idx), results.train_mse(idx), results.test_mse(idx), ratio);
end
fprintf('========================================\n\n');


best_test_mse = min(results.test_mse);
best_test_idx = find(results.test_mse == best_test_mse, 1);


tolerance = 0.05;  
threshold = best_test_mse * (1 + tolerance);


valid_idx = find(results.test_mse <= threshold);
[~, min_M_pos] = min(results.M(valid_idx));
best_idx = valid_idx(min_M_pos);

best_M = results.M(best_idx);
best_params = results.best_params{best_idx};

fprintf('--- Best M Selection (Complexity-Aware) ---\n');
fprintf('Minimum Test MSE: %.6f (at M=%d)\n', best_test_mse, results.M(best_test_idx));
fprintf('Threshold (%.0f%% tolerance): %.6f\n', tolerance*100, threshold);
fprintf('Models within threshold: M = [');
fprintf('%d ', results.M(valid_idx));
fprintf(']\n');
fprintf('Selected M=%d with Test MSE=%.6f\n', best_M, results.test_mse(best_idx));
fprintf('Reason: Lowest complexity among acceptable models\n');
fprintf('==========================================\n\n');



figure('Position', [100 100 1000 400]);
subplot(1,2,1);
plot(results.M, results.train_mse, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(results.M, results.test_mse, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
% Highlight το επιλεγμένο M
plot(best_M, results.test_mse(best_idx), 'g*', 'MarkerSize', 15, 'LineWidth', 2);
grid on; 
legend('MSE Εκπαίδευσης', 'MSE Ελέγχου', ...
       sprintf('Επιλεγμένο μοντέλο (M=%d)', best_M), ...
       'Location', 'best');
xlabel('Αριθμός Gaussian συναρτήσεων (M)');
ylabel('Μέσο Τετραγωνικό Σφάλμα (MSE)');
title('Σφάλμα (MSE) σε συνάρτηση με την πολυπλοκότητα του μοντέλου');



subplot(1,2,2);
x = 1:num_runs;
width = 0.35;
b1 = bar(x - width/2, results.train_mse, width, 'FaceColor', [0.2 0.6 0.8]);
hold on;
b2 = bar(x + width/2, results.test_mse, width, 'FaceColor', [0.8 0.4 0.2]);
set(gca, 'XTick', x, 'XTickLabel', results.M);
xlabel('Αριθμός Gaussian συναρτήσεων (M)');
ylabel('Μέσο Τετραγωνικό Σφάλμα (MSE)');
legend('Train', 'Test', 'Location', 'best');
grid on;
title('Train vs Test MSE');


figure('Position', [100 100 1200 400]);
colors = lines(num_runs);
for idx = 1:num_runs
    plot(results.best_hist{idx}, 'Color', colors(idx,:), 'LineWidth', 1.5);
    hold on;
end
grid on;


legend_labels = cell(num_runs, 1);
for idx = 1:num_runs
    if idx == best_idx
        legend_labels{idx} = sprintf('M=%d (επιλεγμένο)', results.M(idx));
    else
        legend_labels{idx} = sprintf('M=%d', results.M(idx));
    end
end
legend(legend_labels, 'Location', 'best');
xlabel('Γενιά');
ylabel('Καλύτερη Τιμή Συνάρτησης Κόστους (MSE)');
title('Εξέλιξη Γενετικού Αλγορίθμου για διαφορετικό αριθμό Gaussians');

[U1, U2] = meshgrid(linspace(u1_min, u1_max, 50), ...
                    linspace(u2_min, u2_max, 50));
Y_true = f(U1(:), U2(:));
Y_pred = model_predict(U1(:), U2(:), best_params);


zmin = min(Y_true);
zmax = max(Y_true);

figure('Position', [100 100 1200 500]);
subplot(1,2,1);
surf(U1, U2, reshape(Y_true, size(U1)));
title('Πραγματική συνάρτηση');

shading interp;
xlabel('u_1'); ylabel('u_2'); zlabel('y');
zlim([zmin zmax]);  % Κοινό z-range
colorbar;
cb = colorbar;
cb.Label.String = 'Τιμή συνάρτησης y';
view(-37.5, 30);
subplot(1,2,2);
surf(U1, U2, reshape(Y_pred, size(U1)));
title(sprintf('Βέλτιστη προσέγγιση (M=%d, MSE ελέγχου=%.4f)', ...
              best_M, results.test_mse(best_idx)));
shading interp;
xlabel('u_1'); ylabel('u_2'); zlabel('y');
zlim([zmin zmax]);  % Κοινό z-range
colorbar;
cb = colorbar;
cb.Label.String = 'Τιμή συνάρτησης y';
view(-37.5, 30);

fprintf('\n--- Generating detailed plots for best M=%d ---\n', best_M);

yhat_te_best = model_predict(u1_te, u2_te, best_params);


figure('Position', [100 100 600 500]);
scatter(y_te, yhat_te_best, 20, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(y_te) max(y_te)], [min(y_te) max(y_te)], 'r--', 'LineWidth', 2);
grid on;
xlabel('Πραγματική τιμή y (σύνολο ελέγχου)');
ylabel('Προβλεπόμενη τιμή y (σύνολο ελέγχου)');
title(sprintf('Πραγματικές vs προβλεπόμενες τιμές (σύνολο ελέγχου) | M=%d, MSE=%.4f', ...
              best_M, results.test_mse(best_idx)));
legend('Δεδομένα', 'Ιδανική προσαρμογή', 'Location', 'best');
axis equal; axis tight;


figure('Position', [100 100 600 500]);
err_test = y_te - yhat_te_best;
histogram(err_test, 30, 'FaceColor', [0.3 0.7 0.9], 'EdgeColor', 'k');
grid on;
xlabel('Σφάλμα (y - ŷ) στο σύνολο ελέγχου');
ylabel('Συχνότητα');
title('Κατανομή σφάλματος στο σύνολο ελέγχου | M=10');
hold on;
xline(0, 'r--', 'LineWidth', 2);
xline(mean(err_test), 'g--', 'LineWidth', 1.5);
legend('Σφάλμα', 'Μηδενικό σφάλμα', 'Μέση τιμή', 'Location', 'best');

text(mean(err_test), max(ylim)*0.9, ...
     sprintf('Μέση τιμή = %.4f', mean(err_test)), ...
     'HorizontalAlignment', 'center', ...
     'BackgroundColor', 'white');


figure('Position', [100 100 800 500]);
plot(results.best_hist{best_idx}, 'b-', 'LineWidth', 2);
grid on;
xlabel('Γενιά');
ylabel('Καλύτερη τιμή συνάρτησης κόστους (train MSE)');
title('Σύγκλιση Γενετικού Αλγορίθμου | M=10');
max_y = max(results.best_hist{best_idx}(1:min(20, gens)));
ylim([0 max_y*1.1]);

fprintf('\n===========================================\n');
fprintf('All plots generated successfully!\n');
fprintf('Best M = %d with Test MSE = %.6f\n', best_M, results.test_mse(best_idx));
fprintf('===========================================\n');