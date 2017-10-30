c=0
R=zeros(3,3,3)
for lambda=[1,5,10,50,100,200,500,1000,2000,5000,1000]
 [W,TotalTime] = PMFA_L2(saccharomyces.rxnE,S,lambda,10*saccharomyces.L,10*saccharomyces.U,3)
V=varianceCap(saccharomyces.rxnE,W)
M=S*W
Ms=sum(abs(M))
c=c+1
R2(c,1,1)=lambda,
R2(c,2,1:3)=V
R2(c,3,1:3)=Ms
end
c=0
for lambda=[1,5,10,50,100,200,500,1000,2000,5000,10000]
 [W,TotalTime] = PMFA_L1(saccharomyces.rxnE,S,lambda,10*saccharomyces.L,10*saccharomyces.U,3)
V=varianceCap(saccharomyces.rxnE,W)
M=S*W
Ms=sum(abs(M))
c=c+1
R1(c,1,1)=lambda,
R1(c,2,1:3)=V
R1(c,3,1:3)=Ms

end
save('../SuplementaryResult/sachro.mat','R1','R2')