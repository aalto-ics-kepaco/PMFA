# PMFA
Principal Metabolic Flux Mode Analysis
This package contains PMFA and PSMFA code in Matlab. To run this code once need matlab inc=build function optimization routine "quadprog" . IT contains following 2 fucntions

#SPMFA
function [W,runTime]=SPMFA_L2(Einput,S,lambda,Ac,L,U,num,str,ID)
% solve max x^Tcov(E)x - \lambda \|Sx\|_2 
%such that, and L<=x <=U and
% \|x\|_1 <=Ac
%%
%Input:
%Einput= ReactionExpression Matrix : number of reaction x number of samples
%S=Stoichiometric Matrix : number of metabolites x number of reaction
%lambda= model parameter show how strong steadystate constraint will be
% L=lowaer bound of all reaction
% U=upper bound for all reaction
% Ac= sparsity parameter
% str= some string to save intermediate results. Default is "cputime"
% num = how many principal component need to find out. Default 1
% ID = if consider to analysis a subsystem then ID contains lidt pf index
% of target reactions.
%%
%Output:
%W the PMF loadings
% runtime: total time taken

#PMFA
function [W,runTime]=PMFA_L2(Einput,S,lambda,L,U,str,num,ID)
% solve max x^Tcov(E)x - \lambda \|Sx\|_2 such that, and L<=x <=U
%%
%Input:
%Einput= ReactionExpression Matrix : number of reaction x number of samples
%S=Stoichiometric Matrix : number of metabolites x number of reaction
%lambda= model parameter show how strong steadystate constraint will be
% L=lowaer bound of all reaction
% U=upper bound for all reaction
% str= some string to save intermediate results
% num = how many principal component need to find out. Default 1
% ID = if consider to analysis a subsystem then ID contains lidt pf index
% of target reactions.
%%
%Output:
%W the PMF loadings
% runtime: total time taken
