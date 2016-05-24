### $Id: linregEnd.R 53 2007-04-20 12:05:00Z bhm $
# %=============== linregEnd.m ====================
# %  [BetaU,msError,errorObs,Yhat] = linregEnd(Umodel,Y)
# %        performs the part of linregEst that is not performed by linregStart.
# %        See linregEst.m for description.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [BetaU,msError,errorObs,Yhat] = linregEnd(Umodel,Y)
# BetaU = Umodel'*Y;
# Yhat = Umodel*BetaU;
# %%%---%%%  errorObs = Uerror'*Y;
# [U S V] = economySVD(Y-Yhat);
# df_error = size(Umodel,1) - size(Umodel,2);
# errorObs = S*V'; %%%% smartere beregning her???
# if(size(errorObs,1)>df_error)
#     errorObs = errorObs(1:df_error,:);
# end
# msError=sum(errorObs.*errorObs,1)/df_error;
# if(size(errorObs,1)<df_error)
#     errorObs = {errorObs,df_error};  % cell array when incomplete errorObs
# end
#####################################################################
linregEnd = function(Umodel,Y){
BetaU = t(Umodel)%*%Y
Yhat = Umodel%*%BetaU
SVDresid = svd(Y-Yhat)
df_error = nrow(Umodel) - ncol(Umodel)
if(length(SVDresid$d)>1)
   {S = diag(SVDresid$d)}
else
   {S = SVDresid$d}
   #end
errorObs = S%*%t(SVDresid$v) #%%%% smartere beregning her???
if(nrow(errorObs)>df_error)
    errorObs = errorObs[matlabColon(1,df_error),]
#end
msError= colSums(errorObs*errorObs)/df_error;

if(nrow(errorObs)<df_error)
    errorObs = list(errorObs=errorObs,df_error=df_error)  #% cell array when incomplete errorObs
end
list(BetaU=BetaU,msError=msError,errorObs=errorObs,Yhat=Yhat)
}# linregEnd
