permcoefs.plsRbeta <- function(dataset,ind,nt,modele,family=NULL,method="logistic",link="logit",link.phi=NULL,type="ML"){
PLS_beta_wvc(dataY =dataset[,1], dataX=dataset[ind,-1], nt=nt, modele=modele, family=family, keepstd.coeffs=TRUE, method=method, link=link, link.phi=link.phi, type=type)$std.coeffs
}
