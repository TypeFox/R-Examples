### $Id: xy_Obj.R 53 2007-04-20 12:05:00Z bhm $
# %=============== xy_Obj.m ====================
# %  xyObj = xy_Obj(xObj,Y,yNames)
# %      takes an object created by x_Obj as input and
# %      add response values (Y). Further initial computations
# %      for prediction and testing is made.
# %
# %   Output: XyObj is a structure with fields
# %          xObj: same as input
# %        Y(*.*): same as input
# %   yNames{*,1}: same as input
# %     ssTotFull: = sum(sum(Y.^2));
# %         ssTot: = sum(sum(center(Y).^2));
# %           ss: ss's summed over all responses
# %         Beta: regr model: Y = D_om*Beta  (see linregEst)
# %         Yhat: fitted values
# %      YhatStd: stds of fitted values
# %      msError: msError for each response
# %     errorObs: error observations (can be used in multivariate testing)
# %       hypObs: Type II* hypothesis observations (can be used in multivariate testing)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function xyObj = xy_Obj(xObj,Y,yNames)
# xyObj.xObj=xObj;
# xyObj.Y=Y;
# xyObj.yNames=yNames;
# % Continue estimating model where
# %  X = "OM-adjusted model matrix" ,
# %  Y = Y (reponse data)
# [xyObj.Beta,xyObj.msError,xyObj.errorObs,xyObj.Yhat] = linregEnd(xObj.Umodel,Y);
# ss=[];
# xyObj.YhatStd = sqrt(sum(xObj.Umodel.^ 2,2)*xyObj.msError);
# hypObs = cell(size(xObj.D_test));
# for i=1:length(xObj.D_test)
#     hObs = xObj.D_test{i}'*Y;
#     hypObs{i} = hObs;
#     ss = [ss sum(sum(hObs.^2))];
# end
# if(iscell(xyObj.errorObs)) %%%---%%%
#     ss = [ss xyObj.errorObs{2}*sum(xyObj.msError)];
# else
#     ss = [ss size(xyObj.errorObs,1)*sum(xyObj.msError)];
# end
# xyObj.ssTotFull = sum(sum(Y.^2));
# xyObj.ssTot     = sum(sum(center(Y).^2));
# xyObj.ss = ss;
# xyObj.hypObs = hypObs;
##############################################################
xy_Obj = function(xObj,Y){
   xyObj1 = linregEnd(xObj$Umodel,Y)
YhatStd = sqrt( matrix(rowSums(xObj$Umodel^2),,1) %*% xyObj1$msError )
ss=c()
hypObs = vector("list",length(xObj$D_test))
for( i in 1:length(xObj$D_test) ){
  hObs = t(xObj$D_test[[i]])%*%Y
  hypObs[[i]] = hObs
  ss = c(ss,sum(hObs^2))
} #end
if(is.list(xyObj1$errorObs)){
  ss = c(ss, xyObj1$errorObs[[2]]*sum(xyObj1$msError))
}else{
  ss = c(ss,nrow(xyObj1$errorObs)*sum(xyObj1$msError))
}#end
xyObj2 = list(xObj=xObj,Y=Y,YhatStd=YhatStd,hypObs=hypObs,ss=ss,ssTotFull=sum(Y^2),ssTot=sum((stdize(Y, scale = FALSE))^2))
c(xyObj1,xyObj2)
}
