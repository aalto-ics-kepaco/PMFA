function [Edef,Q] = Deflation(E,W)
Edef=E;
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
        %Sdef=Sdef-q*q'*Sdef*q*q';
    % Orthogonal Projection
        Edef= Edef - Edef*q*q'*Edef/(q'*Edef*q);
    % using tech report http://www.eecs.berkeley.edu/Pubs/TechRpts/2012/EECS-2012-99.pdf
    
end
end
