bootplsbeta <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsRbeta, sim="ordinary", stype="i",...){
callplsRbeta <- object$call
#dataset <- cbind(y = eval(callplsRbeta$dataY),eval(callplsRbeta$dataX))
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRbeta$nt)
if(!is.null(callplsRbeta$modele)){modele <- eval(callplsRbeta$modele)} else {modele <- "pls"}
if(!is.null(callplsRbeta$family)){family <- eval(callplsRbeta$family)} else {family <- NULL}
if(!is.null(callplsRbeta$method)){method <- eval(callplsRbeta$method)} else {method <- "logistic"}
if(!is.null(callplsRbeta$link)){link <- eval(callplsRbeta$link)} else {link <- "logit"}
if(!is.null(callplsRbeta$link.phi)){link.phi <- eval(callplsRbeta$link.phi)} else {link.phi <- NULL}
if(!is.null(callplsRbeta$type)){type <- eval(callplsRbeta$type)} else {type <- "ML"}
if(typeboot=="plsmodel"){
temp.bootplsRbeta <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type,...)} else {
boot(data=dataset, statistic=permcoefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi, type=type)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
if(typeboot=="fmodel_np"){
temp.bootplsRbeta <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type,...)} else {
boot(data=dataset, statistic=permcoefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi, type=type)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
if(typeboot=="fmodel_par"){
temp.bootplsRbeta <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type,...)} else {
boot(data=dataset, statistic=permcoefs.plsRbeta, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, method=method, link=link, link.phi=link.phi,type=type)}
indices.temp.bootplsRbeta <- !is.na(temp.bootplsRbeta$t[,1])
temp.bootplsRbeta$t=temp.bootplsRbeta$t[indices.temp.bootplsRbeta,]
temp.bootplsRbeta$R=sum(indices.temp.bootplsRbeta)
temp.bootplsRbeta$call$R<-sum(indices.temp.bootplsRbeta)
return(temp.bootplsRbeta)
}
}