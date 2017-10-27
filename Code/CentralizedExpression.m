function [Ec]=CentralizedExpression(E,axis)
if axis ==1
E=E'
end
for d=1:1:size(E,1)
 Ec(d,:)=E(d,:)-mean(E(d,:));
end
end
