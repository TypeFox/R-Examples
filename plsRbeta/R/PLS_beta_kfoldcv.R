PLS_beta_kfoldcv <- function(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=nrow(dataX), NK=1, grouplist=NULL, random=FALSE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12),weights,method,link=NULL,link.phi=NULL,type="ML") {

    if (missing(weights)) {NoWeights <- TRUE} else {NoWeights <- FALSE}
    res <- NULL
    res$nr <- nrow(dataX)
        if (K > res$nr) {
            cat(paste("K cannot be > than nrow(dataX) =",res$nr,"\n"))
            cat(paste("K is set to", nrow(dataX), "\n"))
            K <- res$nr
            random = FALSE
        }
    call <- match.call(expand.dots=FALSE)
    nt <- eval(nt,parent.frame())
    if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
    if (as.character(call["family"])=="NULL") {
        if (modele=="pls") {call$family<-NULL}
        if (modele=="pls-beta") {family<-NULL}
        if (modele=="pls-glm-Gamma") {call$family<-Gamma(link = "inverse")}
        if (modele=="pls-glm-gaussian") {call$family<-gaussian(link = "identity")}
        if (modele=="pls-glm-inverse.gaussian") {call$family<-inverse.gaussian(link = "1/mu^2")}
        if (modele=="pls-glm-logistic") {call$family<-binomial(link = "logit")}
        if (modele=="pls-glm-poisson") {call$family<-poisson(link = "log")}
        if (modele=="pls-glm-polr") {call$family<-NULL}
    }
    if (!is.null(call$family)) {
        if (is.character(call$family)) {call$family <- get(call$family, mode = "function", envir = parent.frame())}
        if (is.function(call$family)) {call$family <- call$family()}
        if (is.language(call$family)) {call$family <- eval(call$family)}
    }
    if (missing(method)){method<-"logistic"}
    if (is.null(link)){link<-"logit"} else {if(!(link %in% c("logit", "probit", "cloglog", "cauchit", "log", "loglog")) & !is(link,"link-glm")) {link<-"logit"}}
    if (modele=="pls") {cat("\nModel:", modele, "\n\n")}
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {print(family)}
    if (modele %in% c("pls-glm-polr")) {cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}
    if (modele=="pls-beta") {cat("\nModel:", modele, "\n\n");cat("Link:", link, "\n\n");cat("Link.phi:", link.phi, "\n\n");cat("Type:", type, "\n\n")}

    if (as.character(call["tol_Xi"])=="NULL") {call$tol_Xi <- 10^(-12)}
    if (as.character(call["modele"])=="NULL") {call$modele <- "pls"}
    if (as.character(call["limQ2set"])=="NULL") {call$limQ2set <- .0975}
    if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
    folds_kfolds <-vector("list",NK)
    if (NK==1) {respls_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      respls_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          respls_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    if (modele=="pls-beta") {
    if (NK==1) {respls_kfolds_phi <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      respls_kfolds_phi <-vector("list",NK)
        for (jj in 1:NK) {
          respls_kfolds_phi[[jj]] <-vector("list",K)
        }
      }
    }
    }
    if (keepdataY) {
    if (NK==1) {dataY_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      dataY_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          dataY_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    }
    if (keepcoeffs) {
    if (NK==1) {coeffs_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      coeffs_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          coeffs_kfolds[[jj]] <-vector("list",K)
        }
      }
    }
    }
    compl = function (part, set)
    {
        comp = c()
        for (z in set) {
            if (length(which(z == part)) == 0) {
                comp = c(comp, z)
            }
        }
        return(comp)
    }
    for (nnkk in 1:NK) {
            cat(paste("NK:", nnkk, "\n"))
        if (K == res$nr) {
            cat("Leave One Out\n")
            random = FALSE
        }
            cat(paste("Number of groups :", K, "\n"))
        if (!is.list(grouplist)) {
            if (random == TRUE) {
                randsample = sample(1:res$nr, replace = FALSE)
                groups = suppressWarnings(split(randsample, as.factor(1:K)))
            }
            else {
                randsample = sample(1:res$nr, replace = FALSE)
                groups = suppressWarnings(split(randsample, as.factor(1:K)))
                be = 1
                en = 0
                for (z in 1:K) {
                  en = en + length(unlist(groups[z]))
                  groups[z] = list(z = c(be:en))
                  be = en + 1
                }
            }
        }
        else {
            nogroups = grouplist[[nnkk]]
            groups = c()
            for (i in 1:K) groups = c(groups, list(compl(as.vector(unlist(nogroups[i])),
                (1:res$nr))))
        }
        rnames = c()
        for (k in 1:K) rnames = c(rnames, rownames(dataX)[-as.vector(unlist(groups[k]))])
        if (K == 1) {rnames = rownames(dataX)}
        folds = c()
        for (ii in 1:K) {
            nofolds = as.vector(unlist(groups[ii]))
            if (K == 1) {
                folds = c(folds, list(nofolds))
                nofolds = NULL
            }
            else folds = c(folds, list(as.vector(unlist(groups[-ii]))))
            if (K == 1) {
                if(NoWeights){
                temptemp <- PLS_beta_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,link=link,link.phi=link.phi,type=type)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                } else {
                temptemp <- PLS_beta_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,,weights=weights,method=method,link=link,link.phi=link.phi,type=type)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}             
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                }
            else {
                  cat(paste(ii,"\n"))
                  if(NoWeights){
                  temptemp <- PLS_beta_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,link=link,link.phi=link.phi,type=type)
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  } else {
                  temptemp <- PLS_beta_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,weights=weights[-nofolds],method=method,link=link,link.phi=link.phi,type=type) 
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
                  if (keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                  }
        }
        folds_kfolds[[nnkk]]<-folds
    }
results <- list(results_kfolds=respls_kfolds)
if (keepcoeffs) {results$coeffs_kfolds <- coeffs_kfolds}
if (keepfolds) {results$folds <- folds_kfolds}
if (keepdataY) {results$dataY_kfolds <- dataY_kfolds}
if (modele=="pls-beta") {results$results_kfolds_phi <- respls_kfolds_phi}
results$call <- call
results$call$nt <- nt
return(results)
}
