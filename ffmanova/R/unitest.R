### $Id: unitest.R 53 2007-04-20 12:05:00Z bhm $
# %=============== uniTest.m ====================
# % [pValues,stat] = uniTest(modelData,errorData,dfError)
# %     calculates univariate F or t (when DF=1) statistics and
# %     and corresponding p-values.
# %
# %     dfError needed if errordata is incomplete (rows of zeros)
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [pValues,stat] = uniTest(modelData,errorData,dfError)
# dfModel = size(modelData,1);
# if(nargin<3)  %%%-%%% errordata may be incomplete
#     dfError = size(errorData,1);
# end
# if(dfModel==0 | dfError==0)
#    pValues=ones(1,size(modelData,2));
#    stat=zeros(1,size(modelData,2));
#    return;
# end
# errorSS = sum(errorData.^2,1);
# if(dfModel==1) % t-stat
#     stat = modelData ./  sqrt(errorSS/dfError);
#     Fstat = stat.^2;
# else % F-stat
#     modelSS = sum(modelData.^2,1);
#     stat=(dfError/dfModel) * (modelSS ./ errorSS);
#     Fstat = stat;
# end
# %%%  pValues = 1 - cdf('F',Fstat,dfModel,dfError);
# pValues =my_pValueF_(Fstat,dfModel,dfError);
# function pValue = my_pValueF_(f,ny1,ny2)
#      pValue = betainc(ny2*((ny2+ny1*f).^(-1)),ny2/2,ny1/2);
################################################################################
unitest = function(modelData,errorData,dfError=dim(errorData)[1]){
dfModel = dim(modelData)[1];
nYvar  = dim(modelData)[];
if(dfModel==0 | dfError==0){
   pValues = rep(1,nYvar)
   stat = rep(0,nYvar)
   return(list(pValues=pValues,stat=stat))
}#end
errorSS = colSums(errorData^2)
if(dfModel==1){ # t-stat
    stat = modelData /  sqrt(errorSS/dfError);
    Fstat = stat^2;
}else{ # F-stat
    modelSS = colSums(modelData^2)
    stat=(dfError/dfModel) * (modelSS / errorSS);
    Fstat = stat;
 }#end
pValues = pf(Fstat,dfModel,dfError,lower.tail = FALSE)
list(pValues=pValues,stat=stat)
}
