tilt.bootpls <- function(object, typeboot="plsmodel", statistic=coefs.plsR, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, stabvalue=1e6,...){
callplsR <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = eval(callplsR$dataY),eval(callplsR$dataX))
nt <- eval(callplsR$nt)
ifbootfail <- as.matrix(as.numeric(rep(NA, ncol(dataset))))

if(typeboot=="plsmodel"){
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail)}
indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
return(temp.tilt.bootplsR)
}
if(typeboot=="fmodel_np"){
  dataRepYtt <- cbind(y = object$RepY,object$tt)
  wwetoile <- object$wwetoile
  #return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataRepYtt, statistic=coefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataRepYtt, statistic=permcoefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail)}
  indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
  temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
  temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
  temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
return(temp.tilt.bootplsR)
}
if(typeboot=="fmodel_par"){
  temp.tilt.bootplsR <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, ...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail)}
  indices.temp.tilt.bootplsR <- !is.na(temp.tilt.bootplsR$t[,1])
  temp.tilt.bootplsR$t=temp.tilt.bootplsR$t[indices.temp.tilt.bootplsR,]
  temp.tilt.bootplsR$R=sum(indices.temp.tilt.bootplsR)
  temp.tilt.bootplsR$call$R<-sum(indices.temp.tilt.bootplsR)
  return(temp.tilt.bootplsR)
}
}
