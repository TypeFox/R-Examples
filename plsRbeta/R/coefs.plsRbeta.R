coefs.plsRbeta <- function(dataset, ind, nt, modele, family=NULL, method="logistic", link=NULL, link.phi=NULL, type="ML") 
{
    tempcoefs <- try(PLS_beta_wvc(dataY = dataset[ind, 1], dataX = dataset[ind, 
        -1], nt = nt, modele = modele, family=family, method=method, link=link, keepstd.coeffs = TRUE, link.phi=link.phi, type=type)$std.coeffs, silent=TRUE)
    if (is.numeric(tempcoefs)) {
        return(tempcoefs)
    }
    else {
        return(as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset))))))
    }
}
