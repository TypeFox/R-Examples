bootplsglm <- function(object, typeboot="fmodel_np", R=250, statistic=coefs.plsRglmnp, sim="ordinary", stype="i", stabvalue=1e6,verbose=TRUE,...){
callplsRglm <- object$call
maxcoefvalues <- stabvalue*abs(object$Coeffs)
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRglm$nt)
ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset)))))
#dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(!is.null(callplsRglm$family)){family <- eval(callplsRglm$family)} else {family <- NULL}

if(typeboot=="plsmodel"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}

if(typeboot=="fmodel_np"){
dataRepYtt <- cbind(y = object$RepY,object$tt)
wwetoile <- object$wwetoile
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataRepYtt, statistic=coefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail, ...)} else {
boot(data=dataRepYtt, statistic=permcoefs.plsRglmnp, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs)-ncol(object$dataX)))], wwetoile = wwetoile, ifbootfail=ifbootfail)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}

if(typeboot=="fmodel_par"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, maxcoefvalues = maxcoefvalues, ifbootfail=ifbootfail, verbose=verbose)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}
}