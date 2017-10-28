function[Wn]=Normalized(W)
for i=1:1:size(W,2)
    if sqrt(W(:,i)'*W(:,i))>0;
    Wn(:,i)=W(:,i)/sqrt(W(:,i)'*W(:,i));
    else
    Wn(:,i)=W(:,i);
    end
end
end
