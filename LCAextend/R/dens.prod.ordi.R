dens.prod.ordi <-
function(y.x,param,var.list=NULL)
{
    y <- y.x[1:length(param$alpha)]   
    if(length(y.x)==length(param$alpha)) x <- NULL
    else x <- y.x[(length(param$alpha)+1):length(y.x)]
    res <- rep(1,times=nrow(param$alpha[[1]]))
    for(k in 1:nrow(param$alpha[[1]])) for(j in 1:length(param$alpha)) if(!is.na(y[j]))
    {
        S.cov <- length(var.list[[j]])
        S.alp <- ncol(param$alpha[[j]])-S.cov+1

        covar.x <- ifelse(S.cov==0,0,sum(param$alpha[[j]][k,S.alp:(S.alp+S.cov-1)]*x[var.list[[j]]]))
# Version 1.1: correction d'erreurs qui donnaient une probabilité erronée avec 
# des covariables: ajout du champ d'indices 1:(S.alp-1) et passage de covar.x comme
# un argument supplémentaire à p.compute au lieu de l'additionner au vecteur alpha
        res[k] <- res[k]*p.compute(param$alpha[[j]][k,1:(S.alp-1)],covar.x)[y[j]]
    }
    res
}

