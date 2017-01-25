function y=projL1(x,r)

% perform_l1ball_projection - compute the projection on the L1 ball
%
%   y=projL1(x,r);
%
%   x is the projection of y on the set {a \ sum_i |a_i| = r}
%


if size(x,3)>1
    x = x(:,:,1) + 1i*x(:,:,2);
    y = projL1(x,r);
    y = cat(3, real(y), imag(y));
    return;
end
    

if r<0
    error('r should be > 0');
end
if r==0
    y = x*0;
    return;
end

n = length(x(:));
% compute the thresholded L1 norm at each sampled value
s0 = sort( abs(x(:)) );
s = cumsum( s0(end:-1:1) ); s = s(end:-1:1);
s = s - s0 .* (n:-1:1)';
% compute the optimal threshold by interpolation
[i,tmp] = max( find(s>r) );
if isempty(i)
    % it means i=0, s(0)=sum(x)
    ss = sum(x);
    t = ( s(1)-r )/( s(1)-ss ) * (0-s0(1)) + s0(1);
    if t<=0
        y = x; return;
    end
else    
    i = i(end);
    t = ( s(i+1)-r )/( s(i+1)-s(i) ) * (s0(i)-s0(i+1)) + s0(i+1);
end
% do the actual thresholding
y = x*0;
I = find(abs(x)>t);
y(I) = x(I) .* max( 1-t./abs(x(I)), 0 );

return;

%y = x;
%y(abs(x)<t) = 0; 
%y(abs(x)>=t) = y(abs(x)>=t) - sign(x(abs(x)>=t))*t;



end

function g=gradf(K1,cK2,X,a,u)
N=size(K1,1)
g=0*u;
temp =0*(u*u');
for i=1:1:N
    for j=1:1:N
    temp=temp+K1(i,j)*cK2(i,j)*(X(i,:)-X(j,:))'*(X(i,:)-X(j,:));
    end
end
g=-(2*a*u'*temp)';
end
