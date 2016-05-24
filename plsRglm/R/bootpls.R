bootpls <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsR, sim="ordinary", stype="i", stabvalue=1e6, verbose=TRUE,...){
callplsR <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsR$nt)
ifbootfail <- as.matrix(as.numeric(rep(NA, ncol(dataset))))

if(typeboot=="plsmodel"){
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, verbose=verbose, ...)} else {
boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, verbose=verbose)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}

if(typeboot=="fmodel_np"){
dataRepYtt <- cbind(y = object$RepY,object$tt)
wwetoile <- object$wwetoile
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataRepYtt, statistic=coefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail, ...)} else {
boot(data=dataRepYtt, statistic=permcoefs.plsRnp, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, wwetoile = wwetoile, ifbootfail = ifbootfail)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}

if(typeboot=="fmodel_par"){
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, verbose=verbose, ...)} else {
boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, verbose=verbose)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}
}
