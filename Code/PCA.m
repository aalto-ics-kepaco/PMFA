function W=PCA(Einput,num)
% solve max x^TEx 

% centralized
E=CentralizedExpression(Einput);

% covariance
CovE=E*E';

%pca
[W,~]=eigs(CovE,num,'LM');


                   
end
