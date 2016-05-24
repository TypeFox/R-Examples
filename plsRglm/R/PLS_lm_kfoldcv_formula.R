PLS_lm_kfoldcv_formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL,verbose=TRUE) {

    if (missing(data)) {data <- environment(formula)}
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass    
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
#    if (identical(method, "model.frame")) return(mf)
    mt <- attr(mf, "terms")
    attr(mt,"intercept")<-0L
    dataY <- model.response(mf, "any")
    if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
    if (length(dim(dataY)) == 1L) {
        nm <- rownames(dataY)
        dim(dataY) <- NULL
        if (!is.null(nm)) names(dataY) <- nm
        }
    dataX <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
        else matrix(, NROW(dataY), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")

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
    if (keepMclassed) {
    if (NK==1) {Mclassed_kfolds <- list(vector("list", K))}
    else
    {
      if (NK>1)
      {
      Mclassed_kfolds <-vector("list",NK)
        for (jj in 1:NK) {
          Mclassed_kfolds[[jj]] <-vector("list",K)
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
                temptemp <- PLS_lm_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,verbose=verbose)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                } else {
                temptemp <- PLS_lm_wvc(dataY=dataY, dataX=dataX, nt=nt, dataPredictY=dataX,scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,weights=weights,verbose=verbose)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}             
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                if (keepMclassed) {Mclassed_kfolds[[nnkk]][[ii]] = unclass(dataY) !=ifelse(temptemp$valsPredict < 0.5, 0, 1)}
                }
            else {
              if(verbose){cat(paste(ii,"\n"))}
                  if(NoWeights){
                  temptemp <- PLS_lm_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,verbose=verbose)
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  } else {
                  temptemp <- PLS_lm_wvc(dataY=dataY[-nofolds], dataX=dataX[-nofolds,], nt=nt, dataPredictY=dataX[nofolds,], scaleX=scaleX,scaleY=scaleY,keepcoeffs=keepcoeffs,tol_Xi=tol_Xi,weights=weights[-nofolds],verbose=verbose)
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict; attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if(keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
                  if(keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  if(keepMclassed) {Mclassed_kfolds[[nnkk]][[ii]] = unclass(dataY[nofolds]) !=ifelse(temptemp$valsPredict < 0.5, 0, 1)}
                  }
        }
        folds_kfolds[[nnkk]]<-folds
    }
results <- list(results_kfolds=respls_kfolds)
if (keepcoeffs) {results$coeffs_kfolds <- coeffs_kfolds}
if (keepfolds) {results$folds <- folds_kfolds}
if (keepdataY) {results$dataY_kfolds <- dataY_kfolds}
if (keepMclassed) {results$Mclassed_kfolds <- Mclassed_kfolds}
results$call <- call
return(results)
}
