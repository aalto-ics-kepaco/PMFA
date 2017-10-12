function [Ec]=CentralizedExpression(E)
for d=1:1:size(E,1)
 Ec(d,:)=E(d,:)-mean(E(d,:));
end
end