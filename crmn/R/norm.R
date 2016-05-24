##' Fit a model which describes the variation of the labeled internal
##' standards from the biological factors.
##'
##' There is often unwanted variation in among the labeled internal
##' standards which is related to the experimental factors due to
##' overlapping peaks etc. This function fits a model that describes
##' that overlapping variation using a scaled and centered PCA /
##' multiple linear regression model. Scaling is done outside the PCA
##' model.
##' @title Standards model
##' @param object an \code{ExpressionSet} or a \code{matrix}.   Note
##' that if you pass a\code{matrix} have to specify the identity of
##' the standards by   passing the appropriate argument to
##' \code{standards}.
##' @param factors the biological factors described in the pheno data
##' slot if \code{object} is an \code{ExpressionSet} or a design
##' matrix if \code{object} is a \code{matrix}.
##' @param ncomp number of PCA components to use.   Determined by
##' cross-validation if left \code{NULL}
##' @param lg logical indicating that the data should be log
##' transformed 
##' @param fitfunc the function that creates the model fit for
##' normalization, must use the same interfaces as \code{lm}.
##' @param ... passed on to \code{Q2}, \code{pca} (if pcaMethods >
##' 1.26.0), \code{standards} and \code{analytes}
##' @return a list containing the PCA/MLR model, the recommended number of
##' components for that model, the standard deviations and mean values
##' and Q2/R2 for the fit.
##' @seealso \code{makeX}, \code{standardsPred}
##' @examples
##' data(mix)
##' sfit <- standardsFit(mix, "type", ncomp=3)
##' slplot(sfit$fit$pc)
##' ## same thing
##' Y <- exprs(mix)
##' G <- model.matrix(~-1+mix$type)
##' isIS <- fData(mix)$tag == 'IS'
##' sfit <- standardsFit(Y, G, standards=isIS, ncomp=3)
##' @author Henning Redestig
##' @export
standardsFit <- function(object, factors, ncomp=NULL, lg=TRUE, fitfunc=lm, ...) {
  X <- makeX(object, factors)

  if(lg)
    lsta <- log2(t(mexprs(standards(object, ...))))
  else
    lsta <- t(mexprs(standards(object, ...)))
  clsta <- scale(lsta)
  means <- attr(clsta, "scaled:center")
  sds <- attr(clsta, "scaled:scale")
  pfit <- fitfunc(clsta~-1+I(X))
  zbzhate <- cbind(resid(pfit))
  np <- max(1, min(nrow(zbzhate) - 1, ncol(zbzhate) - 1, ncomp))
  hp <- library(help="pcaMethods")$info[[1]]
  ver <- gsub("Version:", "", hp[grep("Version:", hp)])
  if(compareVersion(ver, "1.26.0") == 1)
    pc <- pca(zbzhate, nPcs=np, verbose=FALSE, ...)
  else
    withCallingHandlers(pc <- pca(zbzhate, nPcs=np, method="nipals", verbose=FALSE),
                        warning=pcaMuffle)

  r2 <- NULL -> q2
  best <- min(np, ncomp)
  if(is.null(ncomp) ) {
    withCallingHandlers(q2 <- Q2(pc, zbzhate, nruncv=1, ...), warning=pcaMuffle)
    r2 <- pc@R2
    best <- which.max(q2)
  }
  
  list(fit=list(fit=pfit,pc=pc), ncomp=best, means=means, sds=sds, q2=q2, r2=r2)
}

##' PCA and Q2 issues warnings about biasedness and poorly estimated
##' PCs. The first is non-informative and the poorly estimated PCs will
##' show up as poor overfitting which leads to a choice of fewer PCs
##' i.e. not a problem. This function is mean to muffle those warnings.
##' Only used for version of pcaMethods before 1.26.0.
##' @title Muffle the pca function
##' @param w a warning
##' @return nothing
##' @author Henning Redestig
pcaMuffle <- function(w) if(any(grepl("Precision for components", w),
                                grepl("Validation incomplete", w)))
  invokeRestart( "muffleWarning" )

##' Predicted values for the standards
##' 
##' There is often unwanted variation in among the labeled internal
##' standards which is related to the experimental factors due to
##' overlapping peaks etc. This predicts this effect given a model of
##' the overlapping variance. The prediction is given by
##' \eqn{\hat{X}_{IS}=X_{IS}-X_{IS}B}{XhatIS=XIS-XIS*B}
##' @title Predict effect for new data (or get fitted data)
##' @param model result from \code{standardsFit}
##' @param newdata an \code{ExpressionSet} or \code{matrix} with new
##' data (or the data used to fit the model to get the fitted data)
##' @param factors the biological factors described in the pheno data
##' slot if \code{object} is an \code{ExpressionSet} or a design
##' matrix if \code{object} is a \code{matrix}.
##' @param lg logical indicating that the data should be log
##' transformed 
##' @param ... passed on to \code{standards} and \code{analytes}    
##' @return the corrected data
##' @seealso \code{makeX}, \code{standardsFit}
##' @author Henning Redestig
##' @examples
##' data(mix)
##' fullFit <- standardsFit(mix, "type", ncomp=3)
##' sfit <- standardsFit(mix[,-1], "type", ncomp=3)
##' pred <- standardsPred(sfit, mix[,1], "type")
##' cor(scores(sfit$fit$pc)[1,], scores(fullFit$fit$pc)[1,])
##' ## could just as well have been done as
##' Y <- exprs(mix)
##' G <- model.matrix(~-1+mix$type)
##' isIS <- fData(mix)$tag == 'IS'
##' fullFit <- standardsFit(Y, G, ncomp=3, standards=isIS)
##' sfit    <- standardsFit(Y[,-1], G[-1,], ncomp=3,
##'                         standards=isIS)
##' pred <- standardsPred(sfit, Y[,1,drop=FALSE], G[1,,drop=FALSE], standards=isIS)
##' cor(scores(sfit$fit$pc)[1,], scores(fullFit$fit$pc)[1,])
##' @export
standardsPred <- function(model, newdata, factors, lg=TRUE, ...) {
  X <- makeX(newdata, factors)

  if(model$ncomp == 0)
    return(t(mexprs(standards(newdata, ...))))

  if(lg) 
    lsta <- log2(t(mexprs(standards(newdata, ...))))
  else
    lsta <- t(mexprs(standards(newdata, ...)))

  slsta <- scale(lsta, center=model$means, scale=model$sds)
  ## correct for G
  cslstaE <- slsta - predict(model$fit$fit, data.frame(I(X)))
  ## correct for E, get Tz
  predict(model$fit$pc, cslstaE, pcs=model$ncomp)$scores
}

##' Normalize samples by their weight (as in grams fresh weight)
##'
##' Normalize each sample by dividing by the loaded sample weight. The
##' weight argument is takes from the pheno data (or given as numerical
##' vector with one value per sample). Missing values are not
##' tolerated.
##' @title Normalize by sample weight
##' @param object an \code{ExpressionSet}
##' @param weight a string naming the pheno data column with the
##' weight or a numeric vector with one weight value per sample.
##' @param lg is the assay data already on the log-scale or not. If
##' lg, the weight value is also log-transformed and subtraction is
##' used instead of division.
##' @return the normalized expression set
##' @examples
##' data(mix)
##' w <- runif(ncol(mix),1, 1.3)
##' weightnorm(mix, w)
##' @author Henning Redestig
##' @export
weightnorm <- function(object, weight="weight", lg=FALSE) {
  assd <- exprs(object)
  if(length(weight)[1] & is.character(weight))
    w <- pData(object)[,weight]
  else if (length(weight) == ncol(object) & is.numeric(weight))
    w <- weight
  else
    stop("bad weight")
  if(max(assd, na.rm=TRUE) > 100 & lg)
    warning("found large numbers, is this really on the log-scale?")
  if(lg)
    assd <- sweep(assd, 2, log2(w), "-")
  else
    assd <- sweep(assd, 2, w, "/")
  exprs(object) <- assd
  object
}

##' Fit the parameters for normalization of a metabolomics data set.
##'
##' Normalization is first done by fitting a model and then applying
##' that model either to new data or the same data using
##' \code{normPred}. Five different methods are implemented.
##' \describe{
##' \item{t1}{divide by row-means of the \eqn{L_2}  scaled internal standards}
##' \item{one}{divide by value of a single, user defined, internal standard}
##' \item{totL2}{divide by the square of sums of the full dataset}
##' \item{nomis}{See Sysi-Aho et al.}
##' \item{crmn}{See Redestig et al.}
##' }
##' @title Fit a normalization model
##' @param object an \code{ExpressionSet} or a \code{matrix} (with
##' samples as   columns) in which case the \code{standards} must be
##' passed on via \code{...}
##' @param method chosen normalization method
##' @param one single internal standard to use for normalization
##' @param factors column names in the pheno data slot describing the
##' biological factors. Or a design matrix directly.
##' @param lg logical indicating that the data should be log
##' transformed 
##' @param fitfunc the function that creates the model fit for
##' normalization, must use the same interfaces as \code{lm}. 
##' @param formula if fitfunc has formula interface or not
##' @param ... passed on to \code{standardsFit}, \code{standards},
##'{analytes}
##' @export
##' @return a normalization model
##' @seealso \code{normPred}, \code{standards}, \code{model.matrix}
##' @references Sysi-Aho, M.; Katajamaa, M.; Yetukuri, L. & Oresic,
##' M. Normalization method for metabolomics data using optimal
##' selection of multiple internal standards. BMC Bioinformatics,
##' 2007, 8, 93
##'
##' Redestig, H.; Fukushima, A.; Stenlund, H.; Moritz, T.; Arita, M.;
##' Saito, K. & Kusano, M.  Compensation for systematic
##' cross-contribution improves normalization of mass spectrometry
##' based metabolomics data Anal Chem, 2009, 81, 7974-7980
##' @author Henning Redestig
##' @examples
##' data(mix)
##' nfit <- normFit(mix, "crmn", factors="type", ncomp=3)
##' slplot(sFit(nfit)$fit$pc, scol=as.integer(mix$runorder))
##' ## same thing
##' Y <- exprs(mix)
##' G <- model.matrix(~-1+mix$type)
##' isIS <- fData(mix)$tag == 'IS'
##' nfit <- normFit(Y, "crmn", factors=G, ncomp=3, standards=isIS)
##' slplot(sFit(nfit)$fit$pc, scol=as.integer(mix$runorder))
normFit <- function(object, method, one="Succinate_d4", factors=NULL, lg=TRUE,
                    fitfunc=lm, formula=TRUE,...) {
  if("ccmn" %in% method) 
    method[method == "ccmn"] <- "crmn"
  if(!method %in% c("t1", "avg", "one", "nomis", "totL2", "crmn",
                    "ri", "median", "none"))
    stop("unknown nomrmalization method")

  ana <- t(mexprs(analytes(object, ...)))
  sta <- t(mexprs(standards(object, ...)))

  if(max(ana, na.rm=TRUE) > 100 & !lg)
    warning("found large numbers, should this not be log-transformed?")

  if(lg) {
    lsta <- log2(sta)
    lana <- log2(ana)
  } else {
    lsta <- sta
    lana <- ana
  }    

  wasNa <- !is.finite(lana)
  if((any(!is.finite(lsta)) | any(!is.finite(lsta)) |
      any(!is.finite(lana)) | any(!is.finite(lana))) &
     method %in% c("t1", "nomics", "crmn") & ncol(lsta) > 1) {
    message("Found missing or negative values - applying imputation using ppca")
    lana[!is.finite(lana)] <- NA
    lsta[!is.finite(lsta)] <- NA
    lana <- completeObs(pca(lana, method="ppca"))
    if(all(dim(lsta)) > 1)
      lsta <- completeObs(pca(lsta, method="ppca"))
  }

  means <- colMeans(lana)
  sds <- apply(lana, 2, sd, na.rm=TRUE)

  model <- switch(method,
                  ri = {
                    list(one=one)
                  },
                  one = {
                    list(one=one)
                  },
                  t1 = {
                    stacale <-
                      sqrt(diag(crossprod(sta)) / (nrow(sta)))
                    list(stascale=stacale)
                  },
                  nomis = {
                    ## save the original means
                    mfit <- lm(I(lana)~I(lsta))
                    list(fit=mfit, means=means)
                  },
                  {
                    list(NULL)
                  })
  sfit <- list()
  if(method == "crmn") {
    if(is.null(factors))
      stop("specify biological factors when using crmn")
    sfit <- standardsFit(object, factors, lg=lg, ...)
    if(sfit$ncomp == 0) {
      warning("no structured variance, t1 normalization instead")
      return(Recall(object, "t1", one, factors))
    }
    tz <- standardsPred(sfit, object, factors, lg=lg, ...)
    sclana <- scale(lana)
    if(formula)
      pfit <- fitfunc(sclana~-1+I(tz)) else
    pfit <- fitfunc(y=sclana,x=tz)
    

    model <-  list(fit=pfit,
                  sds=sds, means=means)
  }
  model$wasNa <- wasNa
  return(new("nFit", model=model, sFit=sfit, method=method))
}

##' Normalization methods for metabolomics data
##'
##' Wrapper function for \code{normFit} and \code{normPred}
##' @title Normalize a metabolomics dataset
##' @param object an \code{ExpressionSet}
##' @param method the desired method
##' @param segments normalization in a cross-validation setup, only to use for
##'   validation/QC purposes.
##' @param ... passed on to \code{normFit} and \code{normPred}
##' @author Henning Redestig
##' @seealso \code{normFit}, \code{normPred}
##' @return the normalized dataset
##' @export
##' @examples
##' data(mix)
##' normalize(mix, "crmn", factor="type", ncomp=3)
##' #other methods
##' normalize(mix, "one")
##' normalize(mix, "avg")
##' normalize(mix, "nomis")
##' normalize(mix, "t1")
##' normalize(mix, "ri")
##' normalize(mix, "median")
##' normalize(mix, "totL2")
##' ## can also do normalization with matrices
##' Y <- exprs(mix)
##' G <- with(pData(mix), model.matrix(~-1+type))
##' isIS <- with(fData(mix), tag == "IS")
##' normalize(Y, "crmn", factor=G, ncomp=3, standards=isIS)
normalize <- function(object, method, segments=NULL, ...) {
  if(class(object) == "ExpressionSet")
    pData(object) <- dropunusedlevels(pData(object))

  if(is.null(segments)) {
    fit <- normFit(object, method=method, ...)
    object <- normPred(fit, object, ...)
  }
  else {
    anaobj <- analytes(object, ...)
    exprMat <- mexprs(anaobj)
    for(i in segments) {
      fit <- normFit(object[,-i], method=method, ...)
      exprMat[,i] <- mexprs(normPred(fit, object[,i], ...))
    }
    mexprs(anaobj) <- exprMat
    object <- anaobj
  }
  object
}

##' Predict the normalized data using a previously fitted normalization model.
##'
##' Apply fitted normalization parameters to new data to get normalized data.
##' Current can not only handle matrices as input for methods 'RI' and 'one'.
##' @title Predict for normalization
##' @param normObj the result from \code{normFit}
##' @param newdata an \code{ExpressionSet} or a \code{matrix} (in
##' which case the \code{standards} must be passed on via
##' \code{...}), possibly the same as used to   fit the
##' normalization model in order to get the fitted data.
##' @param factors column names in the pheno data slot describing the
##' biological factors. Or a design matrix.
##' @param lg logical indicating that the data should be log
##' transformed 
##' @param predfunc the function to use to get predicted values from
##' the fitted object (only for crmn)
##' @param ... passed on to \code{standardsPred}, \code{standardsFit},
##'ode{standards}, \code{analytes}
##' @export
##' @return the normalized data
##' @seealso \code{normFit}
##' @examples
##' data(mix)
##' nfit <- normFit(mix, "crmn", factor="type", ncomp=3)
##' normedData <- normPred(nfit, mix, "type")
##' slplot(pca(t(log2(exprs(normedData)))), scol=as.integer(mix$type))
##' ## same thing
##' Y <- exprs(mix)
##' G <- with(pData(mix), model.matrix(~-1+type))
##' isIS <- fData(mix)$tag == 'IS'
##' nfit <- normFit(Y, "crmn", factors=G, ncomp=3, standards=isIS)
##' normedData <- normPred(nfit, Y, G, standards=isIS)
##' slplot(pca(t(log2(normedData))), scol=as.integer(mix$type))
##' @author Henning Redestig
normPred <- function(normObj, newdata, factors=NULL, lg=TRUE, predfunc=predict,...) {
  
  sta <- t(mexprs(standards(newdata, ...)))
  ana <- t(mexprs(analytes(newdata, ...)))
  normedObj <- analytes(newdata, ...)

  normedData <-
    switch(method(normObj),
           t1 = {
             nsta <- sweep(sta, 2, model(normObj)$stascale, "/")
             corrFac <- rowMeans(nsta, na.rm=TRUE)

             sweep(ana, 1, corrFac, "/")
           },
           nomis = {
             if(lg) {
               lsta <- log2(sta)
               lana <- log2(ana)
             }
             else {
               lsta <- sta
               lana <- ana
             }
             ## correct
             lnormedData <-
               lana - predict(model(normObj)$fit, data.frame(I(lsta)))
             ## recenter
             lnormedData <- sweep(lnormedData, 2, model(normObj)$means, "+")

             if(lg)
               exp(lnormedData)
             else
               lnormedData
           },
           crmn = {
             ## filter/log
             if(is.null(factors))
               stop("specify biological factors when using CRMN")
             tz <- standardsPred(sFit(normObj), newdata, factors, lg=lg,...)
             if(lg)
               lana <- log2(ana)
             else
               lana <- ana

             ## center/scale
             lana <-
               scale(lana, center=model(normObj)$means, scale=model(normObj)$sds)
             ## correct
             lnormedData <-
               lana - predfunc(model(normObj)$fit, data.frame(I(tz)))
             ## recenter/rescale
             lnormedData <- sweep(lnormedData, 2, model(normObj)$sds, "*")
             lnormedData <- sweep(lnormedData, 2, model(normObj)$means, "+")
             
             if(lg)
               exp(lnormedData)
             else
               lnormedData
           },
           avg = {
             if(lg) {
               av <- rowMeans(log2(sta), na.rm=TRUE)
               lana <- log2(ana)
               lnormedData <- sweep(lana, 1, av)
               exp(lnormedData)
             }
             else {
               av <- rowMeans(sta, na.rm=TRUE)
               sweep(ana, 1, av)
             }
           },
           one = {
             # TODO this won't work with matrix
             chosen <- match(model(normObj)$one, fData(standards(newdata))$synonym)
             if(is.na(chosen))
               stop(paste("chosen internal standard ", model(normObj)$one,
                          " not found. Choose one of ", 
                          paste(fData(standards(newdata))$synonym, collapse=", ",
                                sep=""),
                          sep=""))
             corrFac <- sta[,chosen]
             if(lg)
               sweep(ana, 1, corrFac, "/")
             else
               sweep(ana, 1, corrFac, "-")
           },
           median = {
             if(lg) {
               normedData <-
                 sweep(ana, 1,
                       apply(ana, 1, median, na.rm=TRUE), "/")
             }
             else{
               normedData <-
                 sweep(ana, 1,
                       apply(ana, 1, median, na.rm=TRUE), "-")
             }
             normedData
             
           },
           totL2 = {
             ana.tmp <- ana
             ana.tmp[is.na(ana.tmp)] <- 0
             cF <- diag(crossprod(t(ana.tmp)))
             if(lg) {
               normedData <-
                 sweep(ana, 1,
                       sqrt(cF / ncol(ana)), "/")
             }
             else {
               normedData <-
                 sweep(ana, 1,
                       sqrt(cF / ncol(ana)), "-")
             }               
             normedData
           },
           ri = {
             if(is.null(as.numeric(fData(analytes(newdata))$RI))) 
               stop("Missing RI data")
             if(all(is.na(as.numeric(fData(standards(newdata))$RI))))
               stop("Missing RI data")
             # TODO this won't work with matrix
             chosen <- grep(model(normObj)$one, fData(standards(newdata))$synonym)
             riMapping <- sapply(as.numeric(fData(analytes(newdata))$RI),
                                 function(x) {
               if(is.na(x))
                 return(chosen)
               which.min((as.numeric(fData(standards(newdata))$RI) - x)^2)
             })
             sapply(1:ncol(ana), function(i) {
               if(lg) 
                 ana[,i] / sta[,riMapping[i]]
               else
                 ana[,i] - sta[,riMapping[i]]
             })
           },
           {
             ana
           })

  ## reset the missing values
  if(any(model(normObj)$wasNa))
    normedData[model(normObj)$wasNa] <- NA

  mexprs(normedObj) <- t(normedData)
  normedObj
}

