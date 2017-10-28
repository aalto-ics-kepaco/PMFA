function W=PCA(Einput,num)
% solve max x^TEx 

% centralized
E=CentralizedExpression(Einput,2);

% covariance
CovE=E*E';

%pca
[W,~]=eigs(CovE,num,'LM');


                   
end
