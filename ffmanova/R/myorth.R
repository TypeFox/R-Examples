### $Id: myorth.R 45 2006-08-15 13:11:29Z bhm $
# %=============== myorth.m ====================
# % U = myorth(X)
# %     Modified version of orth
# %     Uses another value of "tol"
# %     See help rank
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function U = myorth(X)
# tol_ = 1e-9;   %%% !! hard coded constant !! %%%
# if size(X,2)==0
#    U=X;
#    return;
# end
# [U,S,V] = economySVD(X);
# S = diag(S);
# r = length(S);
# meanS = mean(S);
# tol = max(size(U)) * S(1) * tol_; % See help rank
# while S(r) < tol;
#    r=r-1;
# end
# U=U(:,1:r);
# if( size(U,2)==1 & size(X,2)==1 ) % ensure positive correlation
#     if( (X'*U) <0)                % when univariate
#         U = -U;
#     end
# end
#############################################################
myorth = function(X,tol_ = 1e-9){
if(dim(X)[2]==0)
   return(X)
U = svd(X,nv=0)
S = U$d
U = U$u
r = length(S)
# meanS = mean(S); hvorfor var denne med ?????
tol = max(dim(X)) * S[1] * tol_
if(!tol)
   return(X[,numeric(0),drop=FALSE])
while(S[r] < tol)
   r=r-1
U=U[,1:r,drop = FALSE]
if(dim(U)[2] ==1 & dim(X)[2]==1 )
     if( (t(X)%*%U) <0)
        U = -U
U
}# end myorth
