### $Id: linregStart.R 45 2006-08-15 13:11:29Z bhm $
# %=============== linregStart.m ====================
# %  [Umodel,Uerror,VmodelDivS,VextraDivS1] = linregStart(X)
# %        performs the part of linregEst that can be done without knowing Y.
# %        See linregEst.m for description.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [Umodel,VmodelDivS,VextraDivS1] = linregStart(X)
# rank_lim = 1e-9; %%% !! hard coded constant !! %%%
# nObs = size(X,1);
# %%%---%%% [U,S,V] = svd(X);
# [U,S,V] = svd(X,0);
# S = diag(S);
# r = length(S);
# tol = max(size(U)) * S(1) * rank_lim; % See help rank
# while S(r) < tol;
#    r=r-1;
# end
# S=S(1:r);
# Umodel = U(:,1:r);
# %%%---%%%  Uerror = U(:,(r+1):nObs);
# Vmodel = V(:,1:r);
# VmodelDivS = Vmodel ./ (ones(size(V,1),1) * S');
# % above is same as VmodelDivS = Vmodel*inv(diag(S))
# Vextra = V(:,(r+1):size(V,2));
# VextraDivS1 = Vextra/S(1);
#####################################################################
linregStart = function(X,rank_lim = 1e-9){
nObs = nrow(X)
svdX = svd(X,nv=ncol(X))
r = length(svdX$d)
tol = max(dim(X)) * svdX$d[1] * rank_lim
while(svdX$d[r] < tol)
   r=r-1
S = matrix(svdX$d[1:r],1,r) # = S' i matlab
Umodel=svdX$u[,1:r,drop = FALSE]
Vmodel=svdX$v[,1:r,drop = FALSE]
VmodelDivS = Vmodel /  (matrix(1,nrow(Vmodel),1) %*% S)
Vextra =  svdX$v[,matlabColon(r+1,ncol(X)),drop = FALSE]
VextraDivS1 = Vextra/svdX$d[1]
list(Umodel=Umodel,VmodelDivS=VmodelDivS,VextraDivS1=VextraDivS1)
}# end linregStart
