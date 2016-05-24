estim_ncpFAMD<-function (don, ncp.min = 0, ncp.max = 5, method = c("Regularized", 
                                                                    "EM"), method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05, 
                          threshold = 1e-04) 
{
  tab.disjonctif.NA <- function(tab) {
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), 
                                                     "n", "N", "y", "Y"))) 
        dimnames(x) <- list(row.names(tab), paste(nom, 
                                                  levels(moda), sep = "."))
      else dimnames(x) <- list(row.names(tab), levels(moda))
      return(x)
    }
    if (ncol(tab) == 1) 
      res <- modalite.disjonctif(1)
    else {
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
    }
    return(res)
  }
  
  # from missForest package
  prodna<-function (x, noNA){
    n <- nrow(x)
    p <- ncol(x)
    NAloc <- rep(FALSE, n * p)
    NAloc[sample(n * p, floor(n * p * noNA))] <- TRUE
    x[matrix(NAloc, nrow = n, ncol = p)] <- NA
    return(x)
  }
  
#  if(!("numeric"%in%(lapply(don,class)) & "factor"%in%(lapply(don,class)))){stop("Your data set must contain mixed data.")}
  if((sum(sapply(don,is.numeric))==0) || (sum(!sapply(don,is.numeric))==0)){stop("Your data set must contain mixed data.")}
  method <- match.arg(method, c("Regularized", "regularized", 
                                "EM", "em"), several.ok = T)[1]
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", 
                                      "LOO"), several.ok = T)[1]
  method <- tolower(method)
  method.cv <- tolower(method.cv)
  #reagencement des variables
  jeu<-don[,c(which((lapply(don,class))=="numeric"),which((lapply(don,class))=="factor")),drop=F]
  nbquanti<-sum((lapply(don,class))=="numeric")
  
  #suppression niveaux non pris
  jeu <- droplevels(jeu)
  
  vrai.tab = cbind(jeu[,1:nbquanti,drop=F],tab.disjonctif.NA(jeu[,(nbquanti+1):ncol(jeu),drop=F]))
  if (method.cv == "kfold"){
    res = matrix(NA, ncp.max - ncp.min + 1, nbsim)
    pb <- txtProgressBar(
                        min=1/nbsim*100, max=100,style=3)
    for (sim in 1:nbsim){
      jeuNA <- prodna(jeu,pNA)
      
      for (nbaxes in ncp.min:ncp.max) {
        tab.disj.comp <- imputeFAMD(as.data.frame(jeuNA), 
                                   ncp = nbaxes, method = method, threshold = threshold)$tab.disj
        if (sum(is.na(jeuNA)) != sum(is.na(jeu))){ 
          res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - 
                                                   vrai.tab)^2, na.rm = TRUE)/(sum(is.na(tab.disjonctif.NA(jeuNA))) - 
                                                                                 sum(is.na(tab.disjonctif.NA(jeu))))
        }
      }
      setTxtProgressBar(pb, sim/nbsim*100)
    }
    crit = apply(res, 1, mean, na.rm = TRUE)
    names(crit) <- c(ncp.min:ncp.max)
    close(pb)
    result = list(ncp = as.integer(which.min(crit) + ncp.min - 
                                     1), criterion = crit)
    return(result)
  }
  
  if (method.cv == "loo") {
    pb <- txtProgressBar(min = 0, max = 100, style = 3)
    crit <- NULL
    tab.disj.hat <- vrai.tab
    col.in.indicator <- c(0, rep(1,nbquanti),sapply(jeu[,(nbquanti:ncol(jeu)),drop=F], nlevels))
    for (nbaxes in ncp.min:ncp.max) {
      for (i in 1:nrow(jeu)) {
        for (j in 1:ncol(jeu)) {
          if (!is.na(jeu[i, j])) {
            jeuNA <- as.matrix(jeu)
            jeuNA[i, j] <- NA
            tab.disj.hat[i, (cumsum(col.in.indicator)[j] + 
                               1):(cumsum(col.in.indicator)[j + 1])] <- imputeFAMD(as.data.frame(jeuNA), 
                                                                                  ncp = nbaxes, method = method, threshold = threshold)$tab.disj[i, 
                                                                                                                                                 (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 
                                                                                                                                                                                                               1])]
          }
        }
        setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes==(ncp.min:ncp.max))]-1)*nrow(jeu)+i)/(length(ncp.min:ncp.max)*nrow(jeu))*100))    
      }
      crit <- c(crit, mean((tab.disj.hat - vrai.tab)^2, 
                           na.rm = TRUE))
    }
    close(pb)
    names(crit) <- c(ncp.min:ncp.max)
    return(list(ncp = as.integer(which.min(crit) + ncp.min - 
                                   1), criterion = crit))
  }
}