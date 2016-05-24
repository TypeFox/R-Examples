tilt.bootplsglm <- function(object, typeboot="fmodel_np", statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, stabvalue=1e6,...){
callplsRglm <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRglm$nt)
ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset)))))
#dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(!is.null(callplsRglm$family)){family <- eval(callplsRglm$family)} else {family <- NULL}

if(typeboot=="plsmodel"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail,...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
if(typeboot=="fmodel_np"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglmnp} else {permcoefs.plsRglmnp}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  dataRepYtt <- cbind(y = object$RepY,object$tt)
  wwetoile <- object$wwetoile
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataRepYtt, statistic=coefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail, ...)} else {
    tilt.boot(data=dataRepYtt, statistic=permcoefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
if(typeboot=="fmodel_par"){
#return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family))
  temp.tilt.bootplsRglm <- if(!(sim=="permutation")){tilt.boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail,...)} else {
    tilt.boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail)}
  indices.temp.tilt.bootplsRglm <- !is.na(temp.tilt.bootplsRglm$t[,1])
  temp.tilt.bootplsRglm$t=temp.tilt.bootplsRglm$t[indices.temp.tilt.bootplsRglm,]
  temp.tilt.bootplsRglm$R=sum(indices.temp.tilt.bootplsRglm)
  temp.tilt.bootplsRglm$call$R<-sum(indices.temp.tilt.bootplsRglm)
  return(temp.tilt.bootplsRglm)
}
}
