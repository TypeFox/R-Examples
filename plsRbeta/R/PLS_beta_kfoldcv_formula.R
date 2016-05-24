PLS_beta_kfoldcv_formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", family=NULL, K=nrow(dataX), NK=1, grouplist=NULL, random=FALSE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method,control=list(),contrasts=NULL,sparse=FALSE,sparseStop=TRUE,naive=FALSE,link=NULL,link.phi=NULL,type="ML") {

    if (missing(weights)) {NoWeights <- TRUE} else {NoWeights <- FALSE}  
    if (missing(data)) {data <- environment(formula)}
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass    

if (modele %in% c("pls-beta")) {
oformula <- as.formula(formula)
    formula <- Formula::as.Formula(formula)
    if (length(formula)[2L] < 2L) {
        formula <- Formula::as.Formula(formula(formula), ~1)
        simple_formula <- TRUE
    }
    else {
        if (length(formula)[2L] > 2L) {
            formula <- Formula::Formula(formula(formula, rhs = 1:2))
            warning("formula must not have more than two RHS parts")
        }
        simple_formula <- FALSE
    }
    mf$formula <- formula
}

    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {method<-"logistic"}}
}

if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson","pls","pls-glm-polr")) {
    mt <- attr(mf, "terms")
    attr(mt,"intercept")<-0L
    dataY <- model.response(mf, "any")
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
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(dataY)) stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(dataY)), domain = NA)
        }
    }
if (modele %in% "pls-beta") {
mt <- terms(formula, data = data)
mtX <- terms(formula, data = data, rhs = 1L)
mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
attr(mtX,"intercept")<-0L
attr(mtZ,"intercept")<-0L
dataY <- model.response(mf, "any")
dataX <- if (!is.empty.model(mtX)) model.matrix(mtX, mf) else matrix(, NROW(dataY), 0L)
#if (!is.empty.model(mtX)) model.matrix(mtX, mf, contrasts) else matrix(, NROW(dataY), 0L)
dataZ <- if (!is.empty.model(mtZ)) model.matrix(mtZ, mf) else matrix(, NROW(dataY), 0L)
#if (!is.empty.model(mtZ)) model.matrix(mtZ, mf, contrasts) else matrix(, NROW(dataY), 0L)
    if (length(dataY) < 1) 
        stop("empty model")
    if (!(min(dataY) > 0 & max(dataY) < 1)) 
        stop("invalid dependent variable, all observations must be in (0, 1)")
    n <- length(dataY)
    weights <- model.weights(mf)
    if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    if (is.null(weights)) 
        weights <- 1
    if (length(weights) == 1) 
        weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)
    expand_offset <- function(offset) {
        if (is.null(offset)) 
            offset <- 0
        if (length(offset) == 1) 
            offset <- rep.int(offset, n)
        as.vector(offset)
    }
    offsetX <- expand_offset(model.offset(model.part(formula, 
        data = mf, rhs = 1L, terms = TRUE)))
    offsetZ <- expand_offset(model.offset(model.part(formula, 
        data = mf, rhs = 2L, terms = TRUE)))
    if (!is.null(mf$offset)) 
        offsetX <- offsetX + expand_offset(mf[, "(offset)"])
    offset <- list(mean = offsetX, precision = offsetZ)
}
         
    if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
    res <- NULL
    res$nr <- nrow(dataX)
        if (K > res$nr) {
            cat(paste("K cannot be > than nrow(dataX) =",res$nr,"\n"))
            cat(paste("K is set to", nrow(dataX), "\n"))
            K <- res$nr
            random = FALSE
        }
    call <- match.call(expand.dots=FALSE)
    
    
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {method<-"logistic"}}
}

    
    if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
    if (!(modele %in% c("pls","pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson","pls-glm-polr","pls-beta"))) {print(modele);stop("'modele' not recognized")}
    if (!(modele %in% "pls-glm-family") & !is.null(family)) {stop("Set 'modele=pls-glm-family' to use the family option")}
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
if (is.null(link)){link<-"logit"} else {if(!(link %in% c("logit", "probit", "cloglog", "cauchit", "log", "loglog")) & !is(link,"link-glm")) {link<-"logit"}}
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {print(family)}
    if (modele %in% c("pls-glm-polr")) {cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}
    if (modele=="pls-beta") {cat("\nModel:", modele, "\n\n");cat("Link:", link, "\n\n");cat("Link.phi:", link.phi, "\n\n");cat("Type:", type, "\n\n")}
    if (modele=="pls") {cat("\nModel:", modele, "\n\n")}


    if (as.character(call["tol_Xi"])=="NULL") {call$tol_Xi <- 10^(-12)}
    if (as.character(call["modele"])=="NULL") {call$modele <- "pls"}
    if (as.character(call["limQ2set"])=="NULL") {call$limQ2set <- .0975}
    if (as.character(call["start"])=="NULL") {call$start <- NULL}
    if (as.character(call["method"])=="NULL") {call$method <- NULL}
    if (as.character(call["control"])=="NULL") {call$control <- list()}
    if (as.character(call["contrasts"])=="NULL") {call$contrasts <- NULL}
    if (as.character(call["sparse"])=="NULL") {call$sparse <- FALSE}
    if (as.character(call["sparseStop"])=="NULL") {call$sparseStop <- FALSE}
    if (as.character(call["naive"])=="NULL") {call$contrasts <- FALSE}
    if (as.character(call["link"])=="NULL") {call$link <- "logit"}
    if (as.character(call["link.phi"])=="NULL") {call$link.phi <- NULL}
    if (as.character(call["type"])=="NULL") {call$type <- "ML"}
    
    
    
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
                mf2 <- match.call(expand.dots = FALSE)
                m2 <- match(c("nt","modele","family","scaleX","scaleY","keepcoeffs","tol_Xi","weights","subset","start","etastart","mustart","offset","control","method","contrasts","sparse","sparseStop","naive","link","link.phi","type"), names(mf2), 0L)
                mf2 <- mf2[c(1L, m2)]
                mf2[[1L]] <- as.name("PLS_beta_wvc")
                mf2$family <- family
                mf2$weights <- weights
                mf2$dataY <- dataY
                mf2$dataX <- dataX
                mf2$nt <- eval(mf2$nt,parent.frame())
                mf2$dataPredictY <- dataX
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if(match("method",names(call), 0L)==0L){mf2$method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
if(match("method",names(call), 0L)==0L){mf2$method<-"logistic"} else {if(!(call$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {mf2$method<-"logistic"}}
}
                temptemp <- eval(mf2, parent.frame())
                respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                if(!NoWeights) {attr(respls_kfolds[[nnkk]],"XWeights")=weights; attr(respls_kfolds[[nnkk]],"YWeights")=NULL}
                if (keepcoeffs) {coeffskfolds[[nnkk]][[ii]] = temptemp$coeffs}
                if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = NULL}
                if (modele=="pls-beta") {respls_kfolds_phi[[nnkk]][[ii]] = temptemp$valsPredictPhis}
                }
            else {
                  cat(paste(ii,"\n"))
                  mf2 <- match.call(expand.dots = FALSE)
                  m2 <- match(c("nt","modele","family","scaleX","scaleY","keepcoeffs","tol_Xi","weights","subset","start","etastart","mustart","offset","control","method","contrasts","sparse","sparseStop","naive","link","link.phi","type"), names(mf2), 0L)
                  mf2 <- mf2[c(1L, m2)]
                  mf2[[1L]] <- as.name("PLS_beta_wvc")
                  mf2$family <- family
                  mf2$weights <- weights[-nofolds]
                  mf2$dataY <- dataY[-nofolds]
                  mf2$dataX <- dataX[-nofolds,]
                  mf2$nt <- eval(mf2$nt,parent.frame())
                  mf2$dataPredictY <- dataX[nofolds,]
                  temptemp <- eval(mf2, parent.frame())
                  respls_kfolds[[nnkk]][[ii]] <- temptemp$valsPredict
                  if(!NoWeights) {attr(respls_kfolds[[nnkk]][[ii]],"XWeights")=weights[-nofolds]; attr(respls_kfolds[[nnkk]][[ii]],"YWeights")=weights[nofolds]}
                  if (keepcoeffs) {coeffs_kfolds[[nnkk]][[ii]] = temptemp$coeffs}
                  if (keepdataY) {dataY_kfolds[[nnkk]][[ii]] = dataY[nofolds]}
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
results$call$nt <- mf2$nt
return(results)
}
