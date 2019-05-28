function v = varianceCap(E,W)
Ec = CentralizedExpression(E,2);
W = Normalized(W);
v = zeros(1, size(W,2));
for i = 1:size(W,2)
    pc = Ec' * W(:,1:i);
    pc = CentralizedExpression(pc, 1);
    [~,R] = qr(pc);
    if i == 1
        v(i) = sum(R.^2/trace(Ec'*Ec));
    else
        v(i) = sum(diag(R).^2/trace(Ec'*Ec));
    end
end
