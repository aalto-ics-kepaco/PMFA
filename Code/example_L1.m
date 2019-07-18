load("../Data/SaccharomycesFluxomicData.mat");

rand('twister', 2011); % use twister and this seed in order to compare to python

lambda = 2;
n_components = 2;
[W,TotalTime] = PMFA_L1(saccharomyces.rxnE, saccharomyces.StoichiometricMatrix, ...
    lambda, 10*saccharomyces.L, 10*saccharomyces.U, n_components);

%%
V = varianceCap(saccharomyces.rxnE, W);
M = saccharomyces.StoichiometricMatrix * W;
Ms1 = vecnorm(M, 1, 1);
Ms2 = vecnorm(M, 2, 1);
nrm1 = vecnorm(W, 1, 1);
nrm2 = vecnorm(W, 2, 1);

for i = 1:n_components
    fprintf('PMFA_l1: lambda = %g, component %d, cum. exp. var = %.4g, res. l1 = %.4g, res. l2 = %.4g, ||W||_l1 = %.4g, ||W||_l2 = %.4g\n', ...
        lambda, i, V(i), Ms1(i), Ms2(i), nrm1(i), nrm2(i));
end