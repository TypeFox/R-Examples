tilt.bootplsbeta <- function(object, typeboot="plsmodel", statistic=coefs.plsRbeta, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1){
callplsRbeta <- object$call
dataset <- cbind(y = eval(callplsRbeta$dataY),eval(callplsRbeta$dataX))
nt <- eval(callplsRbeta$nt)
if(!is.null(callplsRbeta$modele)){modele <- eval(callplsRbeta$modele)} else {modele <- "pls"}
if(!is.null(callplsRbeta$family)){family <- eval(callplsRbeta$family)} else {family <- NULL}
if(!is.null(callplsRbeta$link)){link <- eval(callplsRbeta$link)} else {link <- "logit"}
if(!is.null(callplsRbeta$link.phi)){link.phi <- eval(callplsRbeta$link.phi)} else {link.phi <- NULL}
if(!is.null(callplsRbeta$type)){type <- eval(callplsRbeta$type)} else {type <- "ML"}
if(typeboot=="plsmodel"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type))
}
if(typeboot=="fmodel_np"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type))
}
if(typeboot=="fmodel_par"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRbeta} else {permcoefs.plsRbeta}, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family, link=link, link.phi=link.phi, type=type))
}
}
