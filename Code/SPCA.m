function W=SPCA(Einput,lambda,num,ID)
% solve max x^TEx - \lambda \|Sx\|_2 
% initialization;

eps=1.0000e-10;
D=size(Einput,1);
N=size(Einput,2);
Nr=length(ID);


% select only mentiond reactions set and centralization of expression

E=zeros(D,N);
E(ID,:)=CentralizedExpression(Einput);

% covariance
CovE=E*E';

% initialization of W
[Ainit,~]=eigs(CovE,num,'LM');
disp('eig complete');



B=0.*Ainit;    


         A=Ainit;
         t=0;
         st=tic,
	diff=norm(A);
         diag(A'*CovE*A)

         while diff>eps
                t=t+1
                A_old=A;
                B_old=B;
                for k=1:1:num;
                    B(:,k)=max(abs(A_old(:,k)'*CovE)-lambda/2,0).* sign(A_old(:,k)'*CovE);
		            end

                [U,~,V]=svd(CovE*B);
                A=U(:,1:num)*V';
                diff=norm(B-B_old)/norm(B+B_old);
                
        W=Normalized(B);
        O.allobj_pca(t,:) = diag(W'*CovE*W);
        O.allobj_sparse(t,:)=sum(abs(W));       
        O.allW(t,:,:)=W;
        O.allobj(t,:)=O.allobj_pca(t,:) - lambda*O.allobj_sparse(t,:);
        disp(['t=',num2str(t),'obj=',num2str(O.allobj(t,:)),', diff=', num2str(diff)]);
	end
        runTime=toc(st);

         W=Normalized(B);

                   
end
