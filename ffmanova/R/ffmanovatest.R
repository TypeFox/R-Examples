### $Id: ffmanovatest.R 53 2007-04-20 12:05:00Z bhm $
# % File: ffmanovatest.m
# %
# % Purpose: Perform Fifty-Fifty MANOVA tests from DFmodel rows of
# %          "model observations" together with DFerror rows of
# %          "error observations".
# %
# % Call from Matlab:
# %     [exVar1,exVar2,dim,dimX,dimY,bufferDim,D,E,A,M,pD,pE,pA,pM] =
# %     ffmanovatest(modelData,errorData,part,partBufDim,minBufDim,
# %            maxBufDim,minErrDf,cp,stand)
# %
# % Input:
# %       modelData(*,*)- The model observations
# %       errorData(*,*)- The error observations
# %       part          - Dimension of Y used in final MANOVA (=dimY)
# %                       test is chosen (when many responses) so that
# %                       part of the variance is explained.
# %                       If part is vector:
# %                       part(1) -> one component.
# %                       part(2) -> two components.
# %                       Last element (= part(j)) -> j,j+1,j+2 ... components
# %       partBufDim    - The part (=dimension) of "free space" used
# %                       as buffer (See computation in code).
# %       minBufDim     - minimum (if possible) dimension of buffer.
# %       maxBufDim     - maximum dimension of buffer.
# %       minErrdf      - minimum dimension of "error space".
# %       cp            - Correction parameter (power) when "few" responses
# %       stand         - Standardization when stand=1
# %                       (e.g. when not same scale for all responses)
# %
# % Output:
# %       exVar1        - Explained Y-variance for dimY components
# %       exVar2        - Explained Y-variance for dimY+bufferdim components
# %       dim           - Dim of final "MANOVA-space"  (= "nObs" - bufferdim)
# %       dimX          - Dim of the X-space in MANOVA (= ordinary DF)
# %       dimY          - Dim of the Y-space in MANOVA (= n components for test)
# %       bufferDim     - Dim of buffer space
# %       D,E,A,M       - Test statistics: D=Wilk , E=Roy , A=Hotell , M=Pillay
# %       pD,pE,pA,pM   - p-values by F-approximation (only lower bound for Roy)
# %
# % NOTE: errorData may be incomplete (rows of zeros)
# %       cell array coding
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [exVar1,exVar2,dim,dimX,dimY,bufferDim,D,E,A,M,pD,pE,pA,pM] = ...
# ffmanovatest(modelData,errorData,part,partBufDim,minBufDim,maxBufDim,minErrDf,cp,stand)
# if(iscell(errorData)) %%%---%%% errorData may not be complete (rows of zeros)
#     dfError = errorData{2};
#     errorData = errorData{1};
#     nZeroRows = dfError-size(errorData,1);
#     %%% errorData = [errorData' zeros(size(errorData,2),nZeroRows)]';
#     min_nZeroRows = size(errorData,2) - size(errorData,1);
#     min_nZeroRows = min(min_nZeroRows,nZeroRows);
#     if(min_nZeroRows>0)
#         errorData = [errorData' zeros(size(errorData,2),min_nZeroRows)]';
#     end
# else
#     dfError = size(errorData,1);
# end
# Y    = [modelData',errorData']';
# dimX = size(modelData,1);
# %  X = zeros(size(Y,1),dimX); %%%% Do not need X explicit %%%%
# %  X(1:dimX,1:dimX) = diag(ones(dimX,1));
# dimFull = dfError + size(modelData,1); %%%---%%%  %dimFull  = size(Y,1);
# dimYfull = min(dimFull,size(Y,2));
# maxDimY  = max(1,min(dimFull-dimX-minErrDf,dimYfull));
# if(dimX==0 | dimYfull==0 | size(errorData,1)==0)
#    exVar1=0;
#    exVar2=0;
#    dim=0;
#    dimY=0;
#    bufferDim=0;
#    D=0;
#    E=0;
#    A=0;
#    M=0;
#    pD=99;
#    pE=99;
#    pA=99;
#    pM=99;
#    return;
# end
# if(stand)
#   for i=1:size(Y,2)
#       Y(:,i) = Y(:,i)/(norm(Y(:,i)));
#   end % standardize
# end
# if(size(Y,2)>size(Y,1))
#     [U,S,V] = economySVD(Y);
# else
#     [U,S,V] = svd(Y); % complete U needed
# end
# if(dimYfull==1 | dimFull==1)
#    ss=S(1,1).^2;
# else
#    ss=diag(S).^2;
# end
# dimY  = 0;
# part_ = 0;
# part_dimY = 1;
# while(part_<part_dimY)
#    dimY = dimY+1;
#    varMod  = sum(ss(1:dimY));
#    varRest = sum(ss) - varMod;
#    if(dimY==maxDimY)
#       part_ = 1000;
#    else
#       factor = sum(((dimY+1):dimFull).^(cp)) /  sum(((dimY+1):dimYfull).^(cp));
#       part_   = varMod / (varMod + factor*varRest);
#    end
#    if(length(part) >= dimY )
#       part_dimY=part(dimY);
#    end
# end
# bufferDim = max(0,floor(1.0001*min(dimYfull-dimY,...
#    partBufDim*(dimFull-dimX-dimY-minErrDf))));
# if(bufferDim>maxBufDim)
#    bufferDim = max(maxBufDim,0);
# end
# if(bufferDim<minBufDim)
#    bufferDim = min(minBufDim,max(0,min(dimYfull-dimY,dimFull-dimX-dimY-minErrDf)));
# end
# exVar1 = varMod/sum(ss);
# exVar2 = sum(ss(1:(dimY+bufferDim)))/sum(ss);
# dim = dimFull - bufferDim;
# XtY = (U(1:dimX,:))';
# XtY = XtY([1:dimY,(dimY+bufferDim+1):end],:); %%%---%%% XtY = XtY([1:dimY,(dimY+bufferDim+1):dimFull],:);
# [XtY,qrR] = qr(XtY,0);
# XtY = XtY(1:dimY,:);
# [U,S,V] = economySVD(XtY);
# if(size(XtY,1)==1 | size(XtY,2)==1)
#    ss=S(1,1).^2;
# else
#    ss=diag(S).^2;
# end
# [D,E,A,M] = multiStatistics(ss);
# [pD,pE,pA,pM] = multiPvalues(D,E,A,M,dim,dimX,dimY);
#####################################################################
# difference from matlab: stand as arg3, + default arguments
ffmanovatest = function(modelData,errorData,stand=0,
                          part=c(0.9,0.5),
                          partBufDim=0.5,
                          minBufDim=0,
                          maxBufDim=100000000,
                          minErrDf=3,
                          cp=-1) {
if(is.list(errorData)) {
    dfError = errorData[[2]];
    errorData = errorData[[1]];
    nZeroRows = dfError- dim(errorData)[1];
    min_nZeroRows = dim(errorData)[2] - dim(errorData)[1];
    min_nZeroRows = min(min_nZeroRows,nZeroRows);
    if(min_nZeroRows>0)
        errorData = rbind(errorData,matrix(0,min_nZeroRows,dim(errorData)[2]))
    #end
} else
    dfError = dim(errorData)[1]
#end
Y  = rbind(modelData,errorData);
dimX = dim(modelData)[1];
dimFull = dfError + dimX;
dimYfull = min(dimFull,dim(Y)[2]);
maxDimY  = max(1,min(dimFull-dimX-minErrDf,dimYfull));
if(dimX==0 | dimYfull==0 | dim(errorData)[1]==0){
   exVar1=0;
   exVar2=0;
   dim=0;
   dimY=0;
   bufferDim=0;
   D=0;
   E=0;
   A=0;
   M=0;
   pD=99;
   pE=99;
   pA=99;
   pM=99;
   return(list(exVar1=exVar1,exVar2=exVar2,dim=dim,dimX=dimX,dimY=dimY,
               bufferDim=bufferDim,
               D=D,E=E,A=A,M=M,pD=pD,pE=pE,pA=pA,pM=pM))
}
#end
if(stand)
  for(i in 1:dim(Y)[2])
      Y[,i] = Y[,i]/(sqrt(sum(Y[,i]^2)))
  #end
#end

Y.svd =svd(Y,nu=dim(Y)[1]);
ss = Y.svd$d^2
dimY  = 0;
part_ = 0;
part_dimY = 1;
while(part_<part_dimY){
   dimY = dimY+1;
   varMod  = sum(ss[1:dimY]);
   varRest = sum(ss) - varMod;
   if(dimY==maxDimY)
      part_ = 1000
   else {
      factor = sum(((dimY+1):dimFull)^(cp)) /  sum(((dimY+1):dimYfull)^(cp));
      part_   = varMod / (varMod + factor*varRest);
   } #end
   if(length(part) >= dimY )
      part_dimY=part[dimY]
   #end
}#end
bufferDim = max(0,floor(1.0001*min(dimYfull-dimY,
   partBufDim*(dimFull-dimX-dimY-minErrDf))));
if(bufferDim>maxBufDim)
   bufferDim = max(maxBufDim,0)
#end
if(bufferDim<minBufDim)
   bufferDim = min(minBufDim,max(0,min(dimYfull-dimY,dimFull-dimX-dimY-minErrDf)))
#end

# Ensure that (dimY+bufferDim) <= rank(Y), where  rank(Y)=myrank(Y,1e-12)
while( bufferDim>0 & (Y.svd$d[dimY+bufferDim]/Y.svd$d[1]) < (max(dim(Y))*1e-12))
  bufferDim = bufferDim-1

exVar1 = varMod/sum(ss);
exVar2 = sum(ss[1:(dimY+bufferDim)])/sum(ss);
dim = dimFull - bufferDim;
XtY = t(Y.svd$u[1:dimX,,drop = FALSE]);
if(bufferDim)
   XtY = XtY[-((dimY+1):(dimY+bufferDim)),,drop = FALSE]
#end
XtY = qr.Q(qr(XtY))
XtY = XtY[1:dimY,,drop = FALSE];
Y.svd =svd(XtY)
ss = Y.svd$d^2
r = multiStatistics(ss);
c(list(exVar1=exVar1,exVar2=exVar2,dim=dim,dimX=dimX,dimY=dimY,bufferDim=bufferDim,D=r$D,E=r$E,A=r$A,M=r$M),
               multiPvalues(r$D,r$E,r$A,r$M,dim,dimX,dimY))
}

# %%%%%%%%%% SUBFUNCTION multiStatistics %%%%%%%%
# function [D,E,A,M] = multiStatistics(ss)
# D = 1;
# A = 0;
# M = 0;
# for i=1:length(ss)
#    D = D*(1-ss(i));
#    if( (1-ss(i))>0 )
#       A = A+(ss(i)/(1-ss(i)));
#    else
#       A = A+ 1/eps;
#    end
#    if(i==1)
#       E = A;
#    end
#    M = M+ss(i);
# end
################################################
multiStatistics = function(ss) {
eps = 2.2204e-016  ### eps in matlab
D = 1;
A = 0;
M = 0;
for( i in 1:length(ss)) {
   D = D*(1-ss[i]);
   if( (1-ss[i])>0 )
      A = A+(ss[i]/(1-ss[i]))
   else
      A = A+ 1/eps;
   #end
   if(i==1)
      E = A;
   #end
   M = M+ss[i];
} #end
list(D=D,E=E,A=A,M=M);
}
# %%%%%%%%%% END SUBFUNCTION multiStatistics %%%%%%%%


# %%%%%%%%%% SUBFUNCTION multiPvalues %%%%%%%%%%
# function [pD,pE,pA,pM] = multiPvalues(D,E,A,M,dim,dimX,dimY)
# p = dimY;
# q = dimX;
# v = dim - dimX;
# s = min(p,q);
# m = (abs(p-q)-1)/2;
# n = (v-p-1)/2;
# r = v - (p-q+1)/2;
# u = (p*q-2)/4;
# t = 1;
# if( (p^2+q^2-5)>0 )
#    t = sqrt( (p^2*q^2-4)/(p^2+q^2-5) );
# end
# r_ = max(p,q);
# if( D^(1/t) <=0 )
#    fD=1e+100;
# else
#    fD = ((1-D^(1/t))/D^(1/t))*((r*t-2*u)/(p*q));
# end
# pD = my_pValueF(fD,p*q,r*t-2*u);
# fE = E*(v-r_+q)/r_;
# pE = my_pValueF(fE,r_,v-r_+q);
# fA = A*2*(s*n+1)/(s^2*(2*m+s+1));
# pA = my_pValueF(fA,s*(2*m+s+1),2*(s*n+1));
# if((s-M)<=0)
#    fM=1e+100;
# else
#    fM = (M/(s-M))*( (2*n+s+1)/(2*m+s+1) );
# end
# pM = my_pValueF(fM,s*(2*m+s+1),s*(2*n+s+1));
#################################################
multiPvalues = function(D,E,A,M,dim,dimX,dimY){
p = dimY;
q = dimX;
v = dim - dimX;
s = min(p,q);
m = (abs(p-q)-1)/2;
n = (v-p-1)/2;
r = v - (p-q+1)/2;
u = (p*q-2)/4;
t = 1;
if( (p^2+q^2-5)>0 )
   t = sqrt( (p^2*q^2-4)/(p^2+q^2-5) );
#end
r_ = max(p,q);
if( D^(1/t) <=0 )
    fD=1e+100
else
   fD = ((1-D^(1/t))/D^(1/t))*((r*t-2*u)/(p*q));
#end
pD = my_pValueF(fD,p*q,r*t-2*u);
fE = E*(v-r_+q)/r_;
pE = my_pValueF(fE,r_,v-r_+q);
fA = A*2*(s*n+1)/(s^2*(2*m+s+1));
pA = my_pValueF(fA,s*(2*m+s+1),2*(s*n+1));
if((s-M)<=0)
   fM=1e+100
else
   fM = (M/(s-M))*( (2*n+s+1)/(2*m+s+1) );
#end
pM = my_pValueF(fM,s*(2*m+s+1),s*(2*n+s+1));
list(pD=pD,pE=pE,pA=pA,pM=pM);
}
# %%%%%%%%%% END SUBFUNCTION multiPvalues %%%%%%%%%%




# %%%%%%%%% SUBFUNCTION my_pValueF %%%%%%%%%%%
# function pValue = my_pValueF(f,ny1,ny2)
# pValue = 100;
# if(isreal(f) & f>0 & ny1>0.9 & ny2>0.9)
#    if(f<1e-13)
#       pValue = 1;
#    else
#       pValue = betainc(ny2/(ny2+ny1*f),ny2/2,ny1/2);
#    end
# end
#################################################
my_pValueF = function(f,ny1,ny2) {
pf(f,ny1,ny2,lower.tail = FALSE)
}
# %%%%%%%%% END SUBFUNCTION my_pValueF %%%%%%%%%%%

