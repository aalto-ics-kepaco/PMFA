function [Covdef,Q] = Deflation(Cov,W)
Covdef=Cov;
Q=[];
T=size(W,2);
N=size(W,1);

for t=1:1:T
    if t==1
        q=W(:,t);
    else
        temp=(eye(N)-Q*Q')*W(:,t);
        if norm(temp)>0.0001
           q=temp/norm(temp);
	else
		q=temp
	end
    end
    Q=[Q,q];
    % Hotelling Deflation 
        %Covdef=Covdef-q*q'*Covdef*q*q';
    % Orthogonal Projection
        Covdef= Covdef - Covdef*q*q'*Covdef/(q'*Covdef*q);
    % using tech report http://www.eecs.berkeley.edu/Pubs/TechRpts/2012/EECS-2012-99.pdf
    
end
end
