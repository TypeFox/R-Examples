PLS_glm_kfoldcv <- function(dataY,dataX,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights, method, verbose=TRUE) {

    if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
    if(!NoWeights){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with weighted PLS\n",sep=""))}} else {NoWeights=TRUE}
    res <- NULL
    res$nr <- nrow(dataX)
        if (K > res$nr) {
            if(verbose){cat(paste("K cannot be > than nrow(dataX) =",res$nr,"\n"))}
            if(verbose){cat(paste("K is set to", nrow(dataX), "\n"))}
            K <- res$nr
            random = FALSE
        }
    call <- match.call(expand.dots=FALSE)
    if(as.character(call["random"])=="NULL"){random=TRUE}
    if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
    if (as.character(call["family"])=="NULL") {
        if (modele=="pls") {call$family<-NULL}
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
    if (modele=="pls") {if(verbose){cat("\nModel:", modele, "\n\n")}}
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {if(verbose){print(family)}}
    if (modele %in% c("pls-glm-polr")) {if(verbose){cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}}

    if (as.character(call["modele"])=="NULL") {call$modele <- "pls"}
    if (as.character(call["limQ2set"])=="NULL") {call$limQ2set <- .0975}
    if (as.character(call["tol_Xi"])=="NULL") {call$tol_Xi <- 10^(-12)}
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
    nnkk=1;while(nnkk<=NK) {
      restartnnkk=FALSE
      if(verbose){cat(paste("NK:", nnkk, "\n"))}
        if (K == res$nr) {
          if(verbose){cat("Leave One Out\n")}
            random = FALSE
        }
      if(verbose){cat(paste("Number of groups :", K, "\n"))}
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
                temptemp <- try(PLS_glm_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele, family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,verbose=verbose),silent=TRUE)
                if(class(temptemp)=="try-error"){restartnnkk=TRUE;break}
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                } else {
                temptemp <- try(PLS_glm_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX, modele=modele, family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,verbose=verbose),silent=TRUE)
                if(class(temptemp)=="try-error"){restartnnkk=TRUE;break}
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}             
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                }
            else {
              if(verbose){cat(paste(ii,"\n"))}
                  if(NoWeights){
                  temptemp <- try(PLS_glm_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,method=method,verbose=verbose),silent=TRUE)
                  if(class(temptemp)=="try-error"){restartnnkk=TRUE;break}
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  } else {
                  temptemp <- try(PLS_glm_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], modele=modele,family=family,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,weights=weights[-nofolds],method=method,verbose=verbose),silent=TRUE) 
                  if(class(temptemp)=="try-error"){restartnnkk=TRUE;break}
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
                  if (keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  }
        }
      if(!restartnnkk){
        folds_kfolds[[nnkk]]<-folds
        nnkk<-nnkk+1
      }
      }
results <- list(results_kfolds=respls_kfolds)
if (keepcoeffs) {results$coeffs_kfolds <- coeffs_kfolds}
if (keepfolds) {results$folds <- folds_kfolds}
if (keepdataY) {results$dataY_kfolds <- dataY_kfolds}
results$call <- call
return(results)
}
