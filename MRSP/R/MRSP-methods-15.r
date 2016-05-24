################################################################################
#####    Methods for MRSP                                                  #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 03.12.2014, 18:09                                  #####
################################################################################

## summary method
summary.MRSP <- function(object, ...){
 zs <- list(DevianceResiduals = summary(residuals(object)),
            BrierScore = object@Brier,
            loglik = object@loglik,
            coef = object@coef,
            edf = object@df,
            AIC = object@AIC,
            BIC = object@BIC)#,
            #iter.count = object@iter.count,
            #best.iter = object@best.iter)
 return(zs)
}

setMethod("summary", "MRSP",
          function(object, ...)
          summary.MRSP(object, ...))


## logLik method
if (!getRversion() >= "2.13.0") {
  setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
} 

setMethod("nobs", signature(object="MRSP"),
function(object, ...){          
  if(all(round(object@weights) == object@weights)){sum(object@weights)}else{nrow(object@y)}
}) 

setMethod("logLik", signature(object="MRSP"),
function(object, ...){
  z <- object@loglik
  attr(z, "df") <- object@df
  attr(z, "nobs") <- nobs(object)
  class(z) <- "logLik"
  z
})

setMethod("AIC", signature(object="MRSP"),
          function(object, ..., k=2)
{
 ll <- logLik(object)
 -2*as.numeric(ll) + k * attr(ll, "df")
})

setMethod("BIC", signature(object="MRSP"),
          function(object, ..., k=2)
{
 ll <- logLik(object)
 -2*as.numeric(ll) + log(nobs(object)) * attr(ll, "df")
})
 

setGeneric("Brier", function(y, mu, weights, ...) standardGeneric("Brier"))
setMethod("Brier", signature(y="matrix"),
function(y, mu, weights = rep(1, nrow(y)), ...){
  weighted.mean(rowMeans((y - mu)^2), weights)
})  
 

## select function: select one MRSP model out of an MRSP.list which is best with
## regards to a certain optimality criterion
## -------- Arguments ------------------------------------------------------- ##
## object:    an object of class "MRSP.list"
## criterion: the criterion based on which the models for different lambdas are
##            compared with each other
## k, type, all the other arguments: arguments to pass to function "cv" if
##            criterion = "cv" is used.
## -------------------------------------------------------------------------- ##
setMethod("select",
          signature(object = "MRSP.list"),
          function(object, criterion = c("AIC", "BIC", "cv"), k = 10,
                   type = "deviance", cvinds = NULL, parallel = TRUE, cores = detectCores(),
                   adaptcv = FALSE, initialseed = sample(1000, size=1), ...)
{
 criterion <- match.arg(criterion)
 
 if(criterion == "cv"){
  cvfit <- cv(object, k=k, type=type, cvinds=cvinds, parallel=parallel, cores=cores,
              adaptcv=adaptcv, initialseed=initialseed, ...)
  if(type == "loglik"){
   bestind <- which.max(cvfit$mean)
  }else{
   bestind <- which.min(cvfit$mean)
  }
 }else if(criterion == "AIC"){
  AICs <- extract(object, "AIC")
  bestind <- which.min(AICs)
 }else if(criterion == "BIC"){
  BICs <- extract(object, "BIC")
  bestind <- which.min(BICs)
 }
 
 out <- object[[bestind]]
 structure(out, "dat" = attr(object, "dat"))
 #return(out)
})


## to be able to use 'update', we need to adjust generic getCall to MRSP:
setMethod("getCall",
          signature(x = "MRSP.list"),
          function(x, ...)
{
 if(!is.null(attr(x, "topcall"))){
  mycall <- attr(x, "topcall")
 }else{
  mycall <- attr(x, "call")
 }
 return(mycall)
})

setMethod("getCall",
          signature(x = "MRSP"),
          function(x, ...)
{
 if(!is.null(attr(x, "topcall"))){
  mycall <- attr(x, "topcall")
 }else{
  mycall <- attr(x, "call")
 }
 return(mycall)
})

## now the update method for MRSP:
setMethod("update",
          signature(object = "MRSP.list"),
          function(object, formula., ..., evaluate = TRUE)
{
 stop("'update' currently does not work for MRSP. This will be fixed in future versions.") 
    if (is.null(call <- getCall(object)))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) {
        oldform <- as.Formula(call$formula)
        #newform <- Formula(formula.)
        call$formula <- as.Formula(Formula:::update.Formula(oldform, formula.)) #update.formula(formula(object), formula.)
    }
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if(!is.null(call$call)) call$call <- NULL
    if (evaluate)
        eval(call, parent.frame())
    else call
})




## the 'show' method 'MRSP.list':
setMethod("show", signature(object = "MRSP.list"),
          function(object)
{
 cat("Model object of class 'MRSP.list'.\n")
 #cat(deparse(substitute(object)), "is a list of", length(object), "objects of class 'MRSP'.\n\n")
 cat("A list of", length(object), "objects of class 'MRSP', with additional attributes 'topcall', 'call' and 'dat'.\n\n")
 if("topcall" %in% names(attributes(object))){
  cat("Topcall:\n")
  print(attr(object, "topcall"))
 }
 cat("\n")
 cat("Model type:\n")
 cat(attr(object, "call")$model@name)
 cat("\n\n")
 cat("Lambda values ranging from ", object[[length(object)]]@lambda, " to ", object[[1]]@lambda, ".\n", sep="")
})


## the 'show' method for 'MRSP':
setMethod("show", signature(object = "MRSP"),
          function(object)
{
 cat("Model object of class 'MRSP'.\n\n")
 if("topcall" %in% names(attributes(object))){
  cat("Call:\n")
  print(attr(object, "topcall"))
 }
 cat("\n")
 cat("Model type:\n")
 cat(object@model@name)
 cat("\n\n")
 cat("Lambda value:\n")
 cat(object@lambda)
 cat("\n\n") 
 cat("LogLik:", object@loglik, "   approx. edf:", object@df, "\n")
 cat("AIC:", object@AIC, "    BIC:", object@BIC, "\n") 
# cat("Type 'R> slotNames([object])' to list all available slots.\n")
})



## a show method for generic 'coef':
setMethod("show", signature(object = "MRSP.coef"),
          function(object)
{ 
 stopifnot(is.list(object))
 hasZ <- length(object) > 1
 haspenindex <- !is.null(attr(object, "penindex"))
 if(haspenindex) penindex <- attr(object, "penindex")
 
 ## get indices of variables with global effect:
 if(haspenindex){
  xglobind <- which(penindex[[1]] %in% c(4, 40, 41, 42))
 }else{
  xglobind <- which(apply(object[[1]], 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5 && !all(u == 0)))
 }
 ## its inverse:
 xcatspecind <- setdiff(seq_len(ncol(object[[1]])), xglobind)
 if(hasZ){
  if(haspenindex){
   zglobind <- which(penindex[[2]] %in% c(2, 20, 21, 22, 23))
  }else{
   zglobind <- which(apply(object[[2]], 2, function(u) diff(range(u)) < .Machine$double.eps ^ 0.5 && !all(u == 0)))
  }
  zcatspecind <- setdiff(seq_len(ncol(object[[2]])), zglobind)
 } 
 ## print output:
 if(length(xcatspecind)){
  cat("Global predictors with category-specific coefficients:\n")
  print(object[[1]][,xcatspecind])
  cat("\n")
 }
 if(length(xglobind)){
  cat("Global predictors with global coefficients:\n")
  print(colMeans(object[[1]][,xglobind]))             ## taking the mean here in case of numerical deviations..
  cat("\n")  
 }
 if(hasZ){
  if(length(zcatspecind)){
   cat("Class-specific predictors with category-specific coefficients:\n")
   print(object[[2]][,zcatspecind])
   cat("\n")
  }
  if(length(zglobind)){
   cat("Category-specific predictors with global coefficients:\n")
   print(colMeans(object[[2]][,zglobind]))             ## taking the mean here in case of numerical deviations..
   cat("\n")  
  }                                                          
 }
})


## the coef method for objects of class MRSP:
setMethod("coef", signature(object = "MRSP"),
          function(object, type=c("original", "standardized", "prethreshold"),
                   simplify=TRUE, ...)
{
 type <- match.arg(type)
 mycoef <-  switch(type,
                   original = object@coef,
                   standardized = object@coef.stand,
                   prethreshold = object@coef.pretres) 

 if(!simplify){
  return(mycoef)
 }else{
  attributes(mycoef) <- list("class" = "MRSP.coef", "penindex" = object@penindex)
  mycoef <- asS4(mycoef)
  show(mycoef)
  invisible(mycoef)
 }
})
coefficients <- coef
 

## a method for fitted values. if convert2hazard=TRUE and a sequential model was
## used, the discrete hazard rates P(Y = r | Y >= r) are returned instead of
## uncoditional probabilities P(Y = r).
setMethod("fitted", signature(object = "MRSP"),
          function(object, convert2hazard = FALSE, ...)
{
 model <- object@model
 if(is.expression(model)) model <- eval(model)
 modelname <- model@name
 
 if(modelname %in% c("Multinomial Logit Model")) out <- mu <- object@mu
 if(modelname %in% c("Sequential Logit Model")){
  if(!convert2hazard){
   out <- mu <- object@mu
  }else{
   mu <- object@mu
   musums <- rowCumsum(mu)
   out <- mu
   K <- ncol(mu)
   out[,seq(2,K)] <- out[,seq(2,K)] / (1-musums[,seq(1,K-1)])
   out
  }
 }
 colnames(out) <- colnames(mu)
 return(out)
})       


## a residuals generic for MRSP:
setMethod("residuals", signature(object = "MRSP"),
          function(object, type = c("deviance", "pearson"), ...)
{
 type <- match.arg(type)
 model <- object@model
 if(is.expression(model)) model <- eval(model)
 modelname <- model@name

 y <- attr(object, "dat")$y
 mu <- fitted(object)
 if(ncol(y) != ncol(mu) | nrow(y) != nrow(mu))
   stop("The dimensions of fitted and original value objects do not match!")

 if(modelname %in% c("Sequential Logit Model")){
  if(!all(rowSums(y) == 1)){
   y <- cbind(y, 1-rowSums(y))
   mu <- cbind(mu, 1-rowSums(mu))
  }
 }

 devvals <- rowSums(y*log(pmax(y, 1e-100)/pmax(mu, 1e-100)))                    ## see "Gerhard Tutz - Regression for Categorical Data (2012, Cambridge University Press)"
 pearsonvals <- rowSums((y - mu)^2/mu)                                          ## for the formulas.

 res <- switch(type,
               "deviance" = sign(rowSums(y-mu))*sqrt(abs(devvals)),
               "pearson"  = sign(rowSums(y-mu))*sqrt(pearsonvals))

 return(res)
})

## predict for MRSP objects. type = response returns probabilities, link returns
## the linear predictors.
setMethod("predict",
          signature(object = "MRSP"),
          function(object, newdata = NULL, type = c("response", "link"),
                   offset, weights, convert2hazard = FALSE, ...)                   ## convert2hazard: if T and model=sequentiallogit() was used, convert probabilities to hazard rates.
{
 type <- match.arg(type)
 if(is.null(newdata)){
  pred <- switch(type, response = object@mu, link = object@eta)
 }else{
  if(is.data.frame(newdata)){
   ## get the new data in the required format:
   mycall <- attr(object, "topcall")
   model <- mycall$model
   if(is.expression(model)) model <- eval(model)
   formula <- mycall$formula  
   K <- ncol(attr(object, "dat")$y)
   if(model@name %in% c("Sequential Logit Model")) K <- K+1
  
   ## fill in some pseudo data for the response since MRSP needs some kind of response variable
   responsename <- formula[[2]]
   dataframe <- get(as.character(mycall$data))
   responseind <- which(names(dataframe) == responsename)
   pseudo.y <- rep(c(diag(K)), (K + nrow(newdata)/K))[seq(1, nrow(newdata))]            # this construction ensures that all categories are occupied.
   newdata <- cbind(pseudo.y, newdata)
   names(newdata)[1] <- as.character(responsename)

   mycall$data <- newdata
   mycall$standardize <- FALSE
   mycall$perform.fit <- FALSE
  
   newmod <- eval(mycall)
   newdata <- newmod$dat[-1]
  }
 ## now we can compute the prediction:
  coef <- object@coef

  if(missing(offset)){ offset = rep(0, nrow(newdata$x)) }
  if(missing(weights)){ weights = rep(1, nrow(newdata$x)) }
  if(ncol(newdata$x) != ncol(coef[[1]])) stop("xnew doesn't match object")
  
  if(is.expression(object@model)) object@model <- eval(object@model)
  invlink <- (object@model)@invlink
  etanew <- updateEta(dat = newdata, coef = coef, weights = weights, offset = offset)
  pred <- switch(type, response = if(!convert2hazard){
                                   invlink(etanew)
                                  }else{
                                   mu <- invlink(etanew)
                                   #onemu <- apply(mu, c(1,2), function(u) 1-u)
                                   #muprods <- rowCumprod(onemu)
                                   #out <- mu
                                   #K <- ncol(mu)
                                   #out[,seq(2,K)] <- out[,seq(2,K)] * muprods[,seq(1,K-1)]
                                   musums <- rowCumsum(mu)
                                   out <- mu
                                   K <- ncol(mu)
                                   out[,seq(2,K)] <- out[,seq(2,K)] / (1-musums[,seq(1,K-1)])
                                   out
                                  }, link = etanew)
 }
 return(pred)
})



## plot functions for objects of class "MRSP.list"
################################################################################
## arguments: 
## x: an object of class "MRSP.list"
## y: a character string giving the name or a number giving the position of the
##    variable for which to plot the paths
## loglambda: if true, the coefficient paths are plotted against the logarithm
##    of 1 + lambda
## stand: if TRUE, the standardized coefficients are plotted.
## lcex: legend cex factor. magnifies (or shrinks) the legend. if set to 0, no
##    legend is plotted.
## smooth: should the resulting paths be smoothed by spline smoothing?
## smoothfactor: factor to multiply the default knot number during spline 
##    smoothing. set this between 0 and 1 to get smoother plots. (the defaults
##    of smooth.spline sometimes dont eliminate all wiggliness from the paths.)
setMethod("plot",
          signature(x = "MRSP.list"),
          function(x, y, type = "l", loglambda = TRUE, stand=FALSE, xlab = NULL, ylab = NULL,
                   col = rainbow(ncol(x[[1]]@dat$y)), main = NULL, lcex = 1,
                   pch = 1, lwd = 2, smooth = TRUE, smoothfactor = 1, lambda = NULL,
                   logbase = exp(1), legendpars = NULL, lambdalwd = NULL, ...)
{
 if(length(y) != 1)
   stop("cannot plot more than one variable at once") 
 lambdas <- extract(x, "lambda")
 lambdas <- sort(lambdas)
 nrlambda <- length(lambdas)
    
 dat <- x[[1]]@dat
 K <- ncol(dat$y)
 hasV <- !is.null(dat$V)
 indg <- x[[1]]@indg
 indcs <- x[[1]]@indcs

 nameslist <- colnames(dat$x)
 if(hasV) nameslist <- c(nameslist, colnames(dat$V[[1]]))
 
 if(!is.numeric(y)){
  what <- pmatch(y, nameslist)
 }else{what <- y}
  
 if(hasV){
  if(what > ncol(dat$x)){
   if((what - ncol(dat$x)) %in% indg){
    coefs <- vector(length = nrlambda) 
    for(i in seq(nrlambda)){
     if(stand){
      coefs[i] <- x[[i]]@coef.stand[[2]][1,(what - ncol(dat$x))]
     }else{
      coefs[i] <- x[[i]]@coef[[2]][1,(what - ncol(dat$x))]
     }
    }
    coefs <- coefs[nrlambda:1]
   }else if((what - ncol(dat$x)) %in% indcs){
    coefs <- matrix(0, nrow = nrlambda, ncol = K)
    for(i in seq(nrlambda)){
     if(stand){
      coefs[i, ] <- x[[i]]@coef.stand[[2]][ ,what - ncol(dat$x)]
     }else{
      coefs[i, ] <- x[[i]]@coef[[2]][ ,what - ncol(dat$x)]
     }
    }
    coefs <- coefs[nrlambda:1,] 
   }
  }      
 }
 if(what <= ncol(dat$x)){
  coefs <- matrix(0, nrow = nrlambda, ncol = K)
  for(i in seq(nrlambda)){
   if(stand){
    coefs[i, ] <- x[[i]]@coef.stand[[1]][ ,what]
   }else{
    coefs[i, ] <- x[[i]]@coef[[1]][ ,what]
   }
  }
  coefs <- coefs[nrlambda:1,]
 }

 if(smooth){
  if(is.vector(coefs)){
   whichnotzero <- which(coefs != 0)
   if(length(unique(whichnotzero)) > 3){
    nk <- ceiling(smoothfactor * n.knots(length(whichnotzero)))
    coefsp <- coefs
    coefsp[whichnotzero] <- smooth.spline(coefs[whichnotzero], nknots = nk)$y
   }else{
    coefsp <- coefs
   }  
  }else{ 
   coefsp <- matrix(0, nrow = nrlambda, ncol = K)
   for(j in seq(K)){
    whichnotzero <- which(coefs[,j] != 0)
    if(length(unique(whichnotzero)) > 3){
     nk <- ceiling(smoothfactor * n.knots(length(whichnotzero)))
     coefsp[whichnotzero,j] <- smooth.spline(coefs[whichnotzero, j], nknots = nk)$y
    }else{
     coefsp[,j] <- coefs[,j]
    }
   }
  }
 }else{coefsp <- coefs}
   
 if(loglambda) lambdas <- log(lambdas + 1, base = logbase)
 if(is.null(main) & !is.numeric(y)){main <- y}
 if(is.null(main) & is.numeric(y)){main <- nameslist[what]}
 #if(is.null(ylab) & (what <= ncol(dat$x))){ylab <- expression(beta)}
 #if(is.null(ylab) & (what > ncol(dat$x))){ylab <- expression(alpha)}
 if(is.null(xlab) & !loglambda){xlab <- expression(lambda)}
 if(is.null(xlab) & loglambda){xlab <- expression(log(1 + lambda ))}
 
 if(!is.array(coefsp)){col <- "black"}

 matplot(lambdas, coefsp, col = col , main = main, type = type, pch = pch,
         lwd = lwd, xlab = xlab, ylab = "", ...)
 if((what <= ncol(dat$x)) & (lcex > 0)){
  if(is.null(colnames(dat$y))){
   legendnames <- seq(2, K + 1)
   for(i in seq(K)){
    legend(x = lambdas[1], y = coefsp[1, i], legend = legendnames[i], bty = "n", xjust = 1, yjust = 0.5, legendpars)
   }
  }else{
   legendnames <- colnames(dat$y)
   if(!is.null(legendpars$x)){lx <- legendpars$x; legendpars$x <- NULL}else{lx <- "topright"}
   if(!is.null(legendpars$lwd)){llwd <- legendpars$lwd; legendpars$lwd <- NULL}else{llwd <- 2*lcex}
   do.call("legend", c(list(x = lx, legend = legendnames, col = col, lwd = llwd, cex = lcex), legendpars))
   #legend(x = lx, legend = legendnames, col = col, lwd = 2*lcex, cex = lcex, legendpars)
  }     
 }
 
 if(is.array(coefsp)){
  identind <- apply(coefsp, 1, function(u){all((abs(u - u[1]) <= 1e-6) & (sign(u) == sign(u[1])))})
  identlist <- list()
  for(i in seq(2, nrlambda)){
   identlist[[i]] <- rep(F, nrlambda)
   if((identind[i] == T) & (identind[i-1] == T)){
    identlist[[i]][c(i-1, i)] <- T
   }
  }
  for(i in seq(2, nrlambda)){
   if(identlist[[i]][i] == T){
    matplot(lambdas[identlist[[i]]], coefsp[identlist[[i]], 1], col = "black", add = T, type = type, pch = pch, lwd = lwd, ...)
   }
  }    
 }
 if(is.array(coefsp) & !(1 %in% x[[1]]@penindex[[1]][what])){
  identind <- apply(coefsp, 1, function(u){any(u == 0)})
  identlist <- list()
  for(i in seq(2, nrlambda)){
   identlist[[i]] <- rep(F, nrlambda)
   if((identind[i] == T) & (identind[i-1] == T)){
   #if(all(identind[seq(i-1, nrlambda)] == T)){
    identlist[[i]][c(i-1, i)] <- T
   }
  }
  for(i in seq(2, nrlambda)){
   if(identlist[[i]][i] == T){
    matplot(lambdas[identlist[[i]]], c(0, 0), col = "black", add = T, type = type, pch = pch, lwd = lwd, ...)
   }
  }    
 } 
 
 if(!is.null(lambda)){
  if(is.null(lambdalwd) & !missing(lwd)) lambdalwd <- lwd
  if(loglambda) abline(v = log(1 + lambda, base = logbase), lty = 2, lwd = lambdalwd)
  else abline(v = lambda, lty = 2, lwd = lambdalwd)
 } 
})  
## todo: legenden automatisch in leeren bereich plotten, eventuell mit der funktion "largest.empty" aus dem Paket Hmisc. 

## Beispiel für Vergrößerung der Grafik:
## normalgroß:
# plot(fit, "Gewerk", loglambda=T, col=partycol, lty=1)
## jetzt alles vergrößert:
# plot(fit, "Gewerk", loglambda=T, col=partycol, lty=1, cex.axis = 1.5, cex.main = 1.5, lcex = 1.5, cex.lab = 1.2, lwd = 2)



##########################################################################################################################
################################################################################
## method for computing bootstrap statistics for MRSP objects
## - object: an object of class "MRSP"
## - dat: data object in the correct form
## - B: number of bootstrap samples
## - bootinds: a prespecified set of bootstrap samples. must be a list of length
##   B.
## - what: the slot of object to sample. if left unspecified, the whole objects
##         are bootstrapped. warning: this can consume extreme amounts of memory!
################################################################################
setMethod("bootstrap",
          signature(object = "MRSP"),
          function(object, dat, what = NULL, B = 50, bootinds = NULL, parallel = TRUE,
                   cores = detectCores(), initialseed = sample(1000, size=1), ...)
{
 objectcall <- object@call
 if(is.expression(objectcall$model)){
  objectcall$model <- eval(objectcall$model)
 } 
 nobs <- nrow(dat$x)
 hasV <- !is.null(dat$V)
 if(missing(what)) what <- NULL
 if(!is.null(what)) what <- as.character(what)
 
 if(missing(bootinds) | is.null(bootinds)){
  set.seed(initialseed)
  bootinds <- list(); length(bootinds) <- B
  for(i in seq(B)){
   sane.ind <- 0
   while(sane.ind < 1){
    bootinds[[i]] <- sample(seq(nobs), replace = T, size = nobs)
    if(length(which(apply(dat$x[bootinds[[i]], ], 2,
                    function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))) == 1){
     if(all(colSums(dat$y[bootinds[[i]], ]) >= 1)){
      sane.ind <- 1
     }
    }
   }
  }
 }else{
  if(length(bootinds) != B) stop("bootinds must be of length B")
 }
 
 # helperfuntions to compute the bootstrap samples:
 bootcorepar <- function(i){
  #source("MRSP-gesamt-8.r")
  library(MRSP)
  bootcore(i = i)
 }
 
 bootcore <- function(i){   
  ind <- bootinds[[i]]

  dat$x <- dat$x[ind,]
  dat$y <- dat$y[ind,]
  if(hasV) dat$V <- lapply(dat$V, function(u){u[ind,]})
  
  objectcall$dat <- quote(dat)
  objectcall$weights <- objectcall$weights[ind]
  objectcall$offset <- objectcall$offset[ind]
  
  m <- eval(objectcall)
  
  if(is.null(what)){
   return(m)
  }else{ 
   return(slot(m, what))
  } 
 }  
 
 ## now the fitting of the bootstrap sample  
 if(parallel){
  cl <- makeCluster(cores)
  clusterExport(cl, ls(envir = sys.frame(sys.nframe())), envir = sys.frame(sys.nframe()))
  bootlist <- parLapply(cl, seq(B), bootcorepar)
  stopCluster(cl)
 }else{
  get(ls(envir = parent.frame()), envir = parent.frame())
  bootlist <- lapply(seq(B), bootcore)
 }
 
 class(bootlist) <- "MRSP.bootstrap"
 return(bootlist)
})  

 
## standard error method for MRSP via bootstrap:
setMethod("se",
          signature(object = "MRSP"),
          function(object, dat, B = 50, bootinds = NULL, method = "bootstrap",
          what = "coef", parallel = TRUE, cores = detectCores(),
          initialseed = sample(1000, size=1), ...)
{
 if(method != "bootstrap") 
   stop(" 'bootstrap' is the only currently supported method for standard error computation")
   
 what <- as.character(what)
 if(!(what == "coef" | what == "coef.stand" | what == "coef.pretres"))
   stop("argument 'what' must either be 'coef', 'coef.stand' or 'coef.pretres'")
   
 bootlist <- bootstrap(object = object, dat = dat, B = B, bootinds = bootinds,
                       what = what, parallel = parallel, cores = cores,
                       initialseed = initialseed, ...)
 
 SE <- slot(object, what)
 for(l in seq_len(length(SE))){
  for(i in seq_len(nrow(SE[[l]]))){
   for(j in seq_len(ncol(SE[[l]]))){
    extractor <- function(u){u[[l]][i,j]}
    SE[[l]][i,j] <- sd(sapply(bootlist, extractor))
   }
  }
 }
 
 if(!is.list(SE)) SE <- list(SE)
 class(SE) <- "MRSP.se"
 return(SE)
})       

## p-value method for MRSP. input: an MRSP object and an object similarly
## structured as object@coef that contains the standard errors.
setMethod("pval",
          signature(object = "MRSP"),
          function(object, SE, what = "coef", ...)
{  
 what <- as.character(what)
 if(!(what == "coef" | what == "coef.stand" | what == "coef.pretres"))
   stop("argument 'what' must either be 'coef', 'coef.stand' or 'coef.pretres'")
   
 coef <- slot(object, what)
 pvals <- Map(function(x,y){x/y}, coef, lapply(SE, function(u){u[u == 0] <- 1;u}))
 pvals <- lapply(pvals, function(u){apply(u, c(1,2), function(v){2*pt(-abs(v),length(object@weights)-object@df)})})
 pvals
})



##############################################################################################

## main function for crossvalidation                         
setMethod("cv",
          signature(object = "MRSP.list"),
          function(object, k = 10, type = "deviance", cvinds = NULL, parallel = TRUE, cores = detectCores(),
                   adaptcv = FALSE, initialseed = sample(1000, size=1), ...)
{
 cores <- min(cores, k)
 if(k < 3) stop("'less than 3'-fold crossvalidation not supported")
 dat <- attr(object[[1]], "dat")
 weights <- object[[1]]@weights
 offset <- object[[1]]@offset
 hasV <- !is.null(dat$V) 
 nobs <- nrow(dat$y)
 P <- ncol(dat$x)
 if(hasV) L <- ncol(dat$V[[1]])
 if(!is.null(attr(object, "call")$control)){      ####
  control <- attr(object, "call")$control
 }else{
  control <- object[[1]]@control
 }
 control@keeparglist <- F
 Refit <- attr(object, "call")$refit
 if(is.null(Refit)) control@keepdat <- F else 
 if(Refit == F) control@keepdat <- F else
 if(Refit == T | Refit == "L") control@keepdat <- T
 control@adaptcv <- F

 if(missing(cvinds) | is.null(cvinds)){
  set.seed(initialseed)
  sane.permutation <- 0
  while(sane.permutation < 1){
   permutation <- sample(seq(nobs))
   boundaries <- numeric(k-1)
   boundaries[1] <- ceiling(nobs / k)
   for(i in seq(2, length(boundaries))){
    if(i%%2 == 1) boundaries[i] <- boundaries[i - 1] + ceiling(nobs / k)
    else boundaries[i] <- boundaries[i - 1] + floor(nobs / k)
   }
   boundaries <- c(0, boundaries)
   boundaries <- c(boundaries, nobs)
   cvinds <- list(); length(cvinds) <- k
   for(i in seq(k)){
    cvinds[[i]] <- permutation[seq((boundaries[i] + 1), boundaries[i+1])]
   }
  
   if(object[[1]]@model@name %in% c("Multinomial Logit Model", "Sequential Logit Model")){ #, "CUB Binomial Logit Model")){
    sane.flag <- numeric(k)
    for(i in seq(k)){
     daty.test <- dat$y[-cvinds[[i]], ]
     if(all(colSums(daty.test) >= 1)){
      if(length(which(apply(dat$x[-cvinds[[i]], ], 2,
                     function(u) diff(range(u)) < .Machine$double.eps ^ 0.5))) <= 1){
       sane.flag[i] <- T
      }
     }
    }
    if(all(sane.flag == T)) sane.permutation <- 1
   }else{
    sane.permutation <- 1
   }
  }
 }else{
  if(length(cvinds) != k) stop("cvinds must be a list of length k")
 }
 
 ## a little helper function to facilitate the use of parLapply:
 cvcorepar <- function(i){
  #source("MRSP-gesamt-8.r")
  library(MRSP)
  cvcore(i = i)
 }
  
 cvcore <- function(i){
  if(!hasV){
   dat <- list(y = dat$y[-cvinds[[i]], ],
               x = dat$x[-cvinds[[i]], , drop = F])
  }else{
   dat <- list(y = dat$y[-cvinds[[i]], ],
               x = dat$x[-cvinds[[i]], , drop = F],
               V = lapply(dat$V, function(u){u[-cvinds[[i]],]}))                
  }
  
  objectcall <- attr(object, "call")  ###attr(object, "call")  ## the last entry of an MRSP.list object was(!) the call
  objectcall$control <- control
  if(is.expression(objectcall$model)){
   objectcall$model <- eval(objectcall$model)
  }
  if(!adaptcv) objectcall$penweights <- object[[1]]@penweights
  objectcall$dat <- quote(dat)
  objectcall$weights <- weights[-cvinds[[i]]]
  objectcall$offset <- offset[-cvinds[[i]]]

  out <- eval(objectcall)
  if(!is.list(out)) out <- list(out)

  out[[1]]@dat <- NULL
  out[[1]]@y <- NULL
  out[[1]]@x.stand <- NULL
  out[[1]]@x.original <- NULL
  if(hasV){
   out[[1]]@V.stand <- NULL
   out[[1]]@V.original <- NULL
  }
  return(out)  
 }

 ## now the fitting of the cv'd models:
 if(parallel){
  cl <- makeCluster(cores)
  #if(sys.parent() > 0) clusterExport(cl, ls(envir = parent.frame()), envir = parent.frame())
  #clusterExport(cl, ls(envir = sys.frame(sys.nframe())), envir = sys.frame(sys.nframe()))
  #clusterExport(cl, list(dat = dat, weights = weights, offset = offset, cvinds = cvinds,
  #                       object = object, hasV = hasV, adaptcv = adaptcv, control = control),
  #              envir = sys.frame(sys.nframe()))
  clusterExport(cl, c("dat", "weights", "offset", "cvinds", as.character(substitute(object)), "hasV", "adaptcv",
                      "control"), envir = sys.frame(sys.nframe()))
  #clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  cvlist <- parLapply(cl, seq(k), cvcorepar)
  stopCluster(cl)
 }else{
  get(ls(envir = parent.frame()), envir = parent.frame())
  cvlist <- lapply(seq(k), cvcore)
 } 

 newdata <- list(); length(newdata) <- k
 for(i in seq(k)){
  if(hasV){
   newdata[[i]] <- list(y = dat$y[cvinds[[i]], ], x = dat$x[cvinds[[i]], ],
                        V = lapply(dat$V, function(u){u[cvinds[[i]],]}),
                        offset = object[[1]]@offset[cvinds[[i]]],
                        weights = object[[1]]@weights[cvinds[[i]]])
  }else{
   newdata[[i]] <- list(y = dat$y[cvinds[[i]], ], x = dat$x[cvinds[[i]], ],
                        offset = object[[1]]@offset[cvinds[[i]]],
                        weights = object[[1]]@weights[cvinds[[i]]])
  }
 }

 loglik <- (object[[1]]@model)@loglik
 invlink <- (object[[1]]@model)@invlink
 nrlambda <- max(1, length(object))
  
 logl <- matrix(nrow = nrlambda, ncol = k)
 dev <- logl
 brierscore <- logl

 for(i in seq(nrow(logl))){
  for(j in seq(ncol(logl))){
   coefij <- cvlist[[j]][[i]]@coef
   etaij <- updateEta(dat = newdata[[j]], coef = coefij,
                      offset = newdata[[j]]$offset, weights = newdata[[j]]$weights)
   muij <- invlink(etaij)
   logl[i, j] <- loglik(y = newdata[[j]]$y, mu = muij, weights = newdata[[j]]$weights)
   dev[i, j] <- 2*(loglik(y = newdata[[j]]$y, mu = newdata[[j]]$y,
                           weights = newdata[[j]]$weights) - logl[i, j])
   brierscore[i, j] <- Brier(y = newdata[[j]]$y, mu = muij, weights = newdata[[j]]$weights)
  } 
 }
 
 means <- switch(type,
                 "loglik" = rowMeans(logl),
                 "deviance" = rowMeans(dev),
                 "brier" = rowMeans(brierscore)
                )
               
 sds <- switch(type,
               "loglik" = apply(logl, 1, sd) / sqrt(k - 1), ## we compute the sd of the estimator for the mean, not the sample sd!!
               "deviance" = apply(dev, 1, sd) / sqrt(k - 1),
               "brier" = apply(brierscore, 1, sd) / sqrt(k - 1)
              )  

 return(list(mean = means, sd = sds, type = type, pred.loglik = logl, pred.deviance = dev, pred.brier = brierscore))
})  


## cv for one MRSP object only instead of a full path:
setMethod("cv",
          signature(object = "MRSP"),
          function(object, k = 5, type = "deviance", cvinds = NULL, parallel = TRUE, cores = detectCores(),
                   adaptcv = FALSE, initialseed = sample(1000, size=1), ...)
{
 cores <- min(cores, k)
 if(k < 3) stop("'less than 3'-fold crossvalidation not supported")


 if(is.null(object@dat) | is.null(attr(object, "dat")))
   stop("in order to cv an MRSP object, you must supply the data in slot 'dat'")
   
 dat <- ifelse(!is.null(object@dat), object@dat, attr(object, "dat"))
 
 objectcall <- object@call
 if(is.expression(objectcall$model)){
  objectcall$model <- eval(objectcall$model)
 }
 object <- list(object)
 objectcall$weights <- NULL                 ## fixme: müsste hier und in cv(object=MRSP.list) eingebaut werden!
 objectcall$offset <- NULL                  ## dito...
 objectcall$coef.init <- objectcall$coef.stand.init <- objectcall$coef.pretres.init <- NULL
 #object[[2]] <- objectcall       ## <- outdated stuff
 object <- structure(object, call = objectcall, dat = dat)
 class(object) <- "MRSP.list"
 object <- asS4(object)
 
 out <- cv(object = object, k = k, type = type, cvinds = cvinds, parallel = parallel, cores = cores,
           adaptcv = adaptcv, initialseed = initialseed, ...)
 return(out)
})


####################################################################################################################################
## the main function for simulations:
################################################################################
## Arguments:
  ## coef: an object of class "MRSP.coef"
  ## calls: a list with the function calls to be used. for example one for group
  ##       lasso and one for the ordinary lasso. the calls must use "dat" for 
  ##       their data argument.
  ## nrep: number of replications used for the simulation
  ## nobs: number of observations used for each replication
  ## generator: a function with argument nobs that creates the design matrix
  ## type: type of measure to be used for model selection in the crossvalidation
  ## predictor.types: list of numeric values that specify the type of predictor
  ##       to be used. 0 means metric predictor, 1 means binary, values > 1 mean
  ##       a categorical predictor with an according number of dummies. note 
  ##       that when using such values, the length of predictor.types might be 
  ##       smaller than the length of the coef object. the connection always is
  ##       "length(coef[[i]]) = sum(predictor.types[[i]])".
  ## k:    how many crossvalidation folds should be used (for each replication)?
  ## refit: should an ML refit based on the active set of the regularized models
  ##       be performed and used in the model selection by cross-validation?
  ## parallel: should parallel computing be employed to use available multicore
  ##       systems?
  ## cores: number of cores for the parallel computations
  ## cvreltol: relative tolerance to be used when selecting the optimal model 
  ##       via crossvalidation. the most sparse model whose mean cv value is not
  ##       more than cvreltol percent worse than the minimum (of the cv means) 
  ##       is chosen as the preferred model. this is due to the fact that the 
  ##       cv-mean function is typically very flat around its min which leads
  ##       to the selection of models which contain significantly more variables
  ##       than necessary while having a cv mean that is something like 0.03-0.2 
  ##       points better than that of much sparser models.
  ## ...:  further arguments to be passed to the functions. this is most likely
  ##       used to pass on arguments for rng, for example mean and sd for rnorm.
  ## -------------------------------------------------------------------------##
  ## output: an object of class MRSP.sim. it is a convoluted list:
  ## the first list corresponds to the different simulation runs
  ## the second one corresponds to the different calls
  ## the third corresponds to the following entries: 1 = bestfit, 2 = pred.err, 
  ##  3 = fit, 4 = cvfit (empty if keepcv = F), 5 = bestm, 6 = truemu (only 
  ##  stored for the first call). 
################################################################################

## function for simulations: 
setMethod("simulation",
          signature(coef = "MRSP.coef"),
function(coef, calls, nrep = 50, nobs = 100, nnewdata = 3*nobs, generator, type = "deviance", k = 5,
         refit = FALSE, keepfit = FALSE, keepcv = TRUE, parallel = TRUE, cores = detectCores(), initialseed = sample(1000, size=1),
         adaptcv = FALSE, cvreltol = 0.01, theta = 0.5,  ...)
{
 if(refit == T) stop("refit not supported yet")
 if(!is.list(calls)) calls <- list(calls)
 
 hasV <- length(coef) > 1
 ncall <- length(calls)
 if(length(cvreltol) == 1 & ncall > 1) cvreltol <- rep(cvreltol, ncall)
 if(length(cvreltol) != ncall) stop("cvreltol has wrong length")
 
 ## a helperfunction for parallelization
 simcorefun <- function(i){
  #source("MRSP-gesamt-8.r")
  library(MRSP)
  
  set.seed(initialseed + i)
  datasim <- generator(nobs)
  model <- calls[[1]]$model
  if(is.expression(model)) model <- eval(model)
  dat <- createResponse(coef = coef, dat = datasim, model = model, theta = theta, ...)
  
  simlist <- list(); length(simlist) <- ncall
  
  for(j in seq(ncall)){
   if(exists("calls[[j]]$refit")) if(calls[[j]]$refit == T) cvreltol[j] <- cvreltol[j] / 4
   fit <- eval(calls[[j]])
   cvfit <- NULL
   if(class(fit) == "MRSP.list"){
    if(type != "AIC" & type != "BIC"){
     cvfit <- eval(call("cv", object = fit, type = type, k = k, parallel = F,
                         adaptcv = adaptcv))
     if(type == "deviance" | type == "brier") bestm <- which.min.cv(cvfit$mean, reltol = cvreltol[j])
     if(type == "loglik") bestm <- which.max.cv(cvfit$mean, reltol = cvreltol[j])
    } 
    if(type == "BIC") bestm <- which.min.cv(extract(fit, "BIC"), reltol = cvreltol[j])
    if(type == "AIC") bestm <- which.min.cv(extract(fit, "AIC"), reltol = cvreltol[j]) 
    bestfit <- fit[[bestm]]
   }else{
    cvfit <- fit
    bestm <- 1
    bestfit <- fit
   }
   
   if(!keepfit) fit <- NULL
   if(!keepcv) cvfit <- NULL

   if(j == 1){
    trueeta <- updateEta(dat = dat, coef = coef, offset = bestfit@offset, weights = bestfit@weights)
    truemu <- (bestfit@model)@invlink(trueeta)
   }else{
    truemu <- NULL
   }
   
   ## prediction error
   loglik <- (bestfit@model)@loglik
   newdatsim <- generator(nnewdata)
   newdat <- createResponse(coef = coef, dat = newdatsim, ...)
   
   pred.eta <- updateEta(dat = newdat, coef = bestfit@coef, offset = rep(0, nrow(newdat$y)), weights = rep(1, nrow(newdat$y)))
   pred.mu <- (bestfit@model)@invlink(pred.eta)
   
   pred.err <- switch(type,
                      "loglik" = loglik(y = newdat$y, mu = pred.mu, weights = rep(1, nrow(newdat$y))),
                      "deviance" = 2*(loglik(y = newdat$y, mu = newdat$y, weights = rep(1, nrow(newdat$y)))
                                      - loglik(y = newdat$y, mu = pred.mu, weights = rep(1, nrow(newdat$y)))),
                      "brier" = Brier(y = newdat$y, mu = pred.mu, weights = rep(1, nrow(newdat$y))),
                      "BIC" = 2*(loglik(y = newdat$y, mu = newdat$y, weights = rep(1, nrow(newdat$y)))
                                      - loglik(y = newdat$y, mu = pred.mu, weights = rep(1, nrow(newdat$y)))), ## same as deviance
                      "AIC" = 2*(loglik(y = newdat$y, mu = newdat$y, weights = rep(1, nrow(newdat$y)))
                                      - loglik(y = newdat$y, mu = pred.mu, weights = rep(1, nrow(newdat$y)))) ## same as deviance  
                     )   
     
   simlist[[j]] <- list(bestfit = bestfit, pred.err = pred.err, fit = fit, cvfit = cvfit, bestm = bestm, truemu = truemu)
  }
  
  return(simlist)
 }  
 
 ## now the fitting of the models for each simulation replication:
 if(parallel){
  cl <- makeCluster(cores)
  clusterExport(cl, ls(envir = sys.frame(sys.nframe())), envir = sys.frame(sys.nframe()))
  #clusterExport(cl, ls(envir = .GlobalEnv), envir = .GlobalEnv)
  simresult <- parLapply(cl, seq(nrep), simcorefun)
  stopCluster(cl)
 }else{
  simresult <- lapply(seq(nrep), simcorefun)
 } 
 
 simresult[[nrep + 1]] <- coef
 
 class(simresult) <- "MRSP.sim"
 return(simresult)
})           


################################################################################

## refitting method for MRSP objects
setMethod("refit",
          signature(object = "MRSP", arglist = "missing"),
function(object, arglist, ...)
{
 objectcall <- object@call
 if(is.expression(objectcall$model)){
  objectcall$model <- eval(objectcall$model)
 }
 objectcall$refit <- F
 refitold <- object@refit
 ridgestabilold <- (object@control)@ridgestabil
 ridgestabilrfold <- (object@control)@ridgestabilrf
 
 objectcall$coef.init <- object@coef
 objectcall$coef.stand.init <- object@coef.stand
 objectcall$coef.pretres.init <- object@coef.pretres
 
 guessed.active <- object@guessed.active
 guessed.active.coef <- object@guessed.active.coef
 guessed.active.groupdiff <- object@guessed.active.groupdiff
 guessed.active.diff <- object@guessed.active.diff
 penweights <- penweightsold <- object@penweights
 
 penweights[[1]] <- lapply(guessed.active, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[2]] <- lapply(guessed.active.coef, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[3]] <- lapply(guessed.active.groupdiff, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[4]] <- lapply(guessed.active.diff, function(u) ifelse(u == T | u == 1, 0, 1e28))
 #sl.indg <- ##
 
 
 objectcall$penweights <- penweights
 
 out <- eval(objectcall, parent.frame(2))
 out@penweights <- penweightsold
 out@call$penweights <- penweightsold
 out@call$refit <- refitold
 out@call$control@ridgestabil <- ridgestabilold
 
 return(out)
}) 


## a helper refitting method to be used in refit-MRSP.list
setMethod("refit",
          signature(object = "MRSP", arglist = "list"),
function(object, arglist, ...)
{
 objectcall <- object@call
 if(is.expression(objectcall$model)){
  objectcall$model <- eval(objectcall$model)
 }
 arglist$refit <- F
 refitold <- object@refit
 ridgestabilold <- (object@control)@ridgestabil
 ridgestabilrfold <- (object@control)@ridgestabilrf
 
 guessed.active <- object@guessed.active
 guessed.active.coef <- object@guessed.active.coef
 guessed.active.groupdiff <- object@guessed.active.groupdiff
 guessed.active.diff <- object@guessed.active.diff
 penweights <- penweightsold <- object@penweights
 
 penweights[[1]] <- lapply(guessed.active, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[2]] <- lapply(guessed.active.coef, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[3]] <- lapply(guessed.active.groupdiff, function(u) ifelse(u == T | u == 1, 0, 1e28))
 penweights[[4]] <- lapply(guessed.active.diff, function(u) ifelse(u == T | u == 1, 0, 1e28))
 
 objectcall$penweights <- arglist$penweights <- penweights

 if(ridgestabilrfold == T) arglist$control@ridgestabil <- T

 ## only keep the elements of arglist that are used in objectcall, i.e. the call
 ## that created object.
 arglist <- arglist[charmatch(intersect(names(objectcall), names(arglist)), names(arglist))]
 
 object <- eval(as.call(c(as.symbol("MRSP.fit"), arglist)))
 object@penweights <- penweightsold
 object@call$penweights <- penweightsold
 object@call$refit <- refitold
 object@call$control@ridgestabil <- ridgestabilold
 
 return(object)
}) 


## refitting method for MRSP.list objects
setMethod("refit",                                                                            
          signature(object = "MRSP.list", arglist = "missing"),
function(object, arglist, ...)                                                           
{
 #if(is.call(attr(object, "call"))) object <- object[-length(object)]       ## since the calls are now stored in an attribute "call", this line is no longer needed

 argl1 <- object[[1]]@arglist
 object[[1]]@arglist <- NULL
 argl1$coef.init <- argl1$coef
 argl1$coef.stand.init <- argl1$coef.stand
 argl1$coef.pretres.init <- argl1$coef.pretres
 
 object[[1]] <- refit(object[[1]], argl1)#object[[1]]
 
 if(length(object) > 1){
  for(pos in seq(2, length(object))){
   argl <- object[[pos]]@arglist
   object[[pos]]@arglist <- NULL
   argl$dat <- object[[1]]@dat
   argl$coef <- argl$coef.init <- object[[pos-1]]@coef
   argl$coef.stand <- argl$coef.stand.init <- object[[pos-1]]@coef.stand
   argl$coef.pretres <- argl$coef.pretres.init <- object[[pos-1]]@coef.pretres

   object[[pos]] <- refit(object[[pos]], argl)
  }
 }
 
 if(length(object) == 1) object <- object[[1]]
 if(length(object) > 1){
  class(object) <- "MRSP.list"
 }
 return(object)
}) 