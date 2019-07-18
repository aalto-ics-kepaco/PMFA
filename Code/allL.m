load("../Data/SaccharomycesFluxomicData.mat");
S = saccharomyces.StoichiometricMatrix;

c = 0;
R2 = zeros(3,3,3);
for lambda = [1,5,10,50,100,200,500,1000,2000,5000,10000]
    [W, TotalTime] = PMFA_L2(saccharomyces.rxnE, saccharomyces.StoichiometricMatrix, ...
        lambda, 10*saccharomyces.L, 10*saccharomyces.U, 3);
    V = varianceCap(saccharomyces.rxnE, W);
    M = S * W;
    Ms = sum(abs(M)); % the residual of S*v for each EFM
    c = c+1;
    R2(c,1,1) = lambda;
    R2(c,2,1:3) = V;
    R2(c,3,1:3) = Ms;
end
%%
explained_variance = R2(:,2,3);
residual_nullspace = sum(R2(:, 2, :), 3);
scatter(explained_variance, residual_nullspace);
for i = 1:size(R2, 1)
    text(explained_variance(i), residual_nullspace(i), num2str(R2(i, 1, 1)));
end

%%
c = 0;
R1 = zeros(3,3,3);
for lambda = [1,5,10,50,100,200,500,1000,2000,5000,10000]
    [W, TotalTime] = PMFA_L1(saccharomyces.rxnE, saccharomyces.StoichiometricMatrix, ...
        lambda, 10*saccharomyces.L, 10*saccharomyces.U, 3);
    V = varianceCap(saccharomyces.rxnE, W);
    M = S * W;
    Ms = sum(abs(M));
    c = c+1;
    R1(c,1,1) = lambda;
    R1(c,2,1:3) = V;
    R1(c,3,1:3) = Ms;
end

%%
explained_variance = R1(:,2,3);
residual_nullspace = sum(R1(:, 2, :), 3);
scatter(explained_variance, residual_nullspace);
for i = 1:size(R1, 1)
    text(explained_variance(i), residual_nullspace(i), num2str(R1(i, 1, 1)));
end

%%
save('../SuplementaryResult/sachro.mat', 'R1', 'R2')