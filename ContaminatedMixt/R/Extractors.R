.ICnames <- function(x=NULL){
  v <- c("AIC", "AICc", "AICu", "AIC3", "AWE","BIC", "CAIC", "ICL")
  if (is.null(x)) return(v)
  else return (v[x])
}
.ModelNames <- function (model){
  type <- switch(EXPR = as.character(model), E = "univariate, equal variance", 
                 V = "univariate, unequal variance", EII = "spherical, equal volume", 
                 VII = "spherical, varying volume", EEI = "diagonal, equal volume and shape", 
                 VEI = "diagonal, equal shape", EVI = "diagonal, equal volume, varying shape", 
                 VVI = "diagonal, varying volume and shape", EEE = "ellipsoidal, equal volume, shape and orientation", 
                 VEE = "ellipsoidal, equal shape and orientation", EVE = "ellipsoidal, equal volume and orientation", 
                 VVE = "ellipsoidal, equal orientation", EEV = "ellipsoidal, equal volume and shape", 
                 VEV = "ellipsoidal, equal shape", EVV = "ellipsoidal, equal volume", 
                 VVV = "ellipsoidal, varying volume, shape, and orientation", 
                 X = "univariate normal", XII = "spherical multivariate normal", 
                 XXI = "diagonal multivariate normal", XXX = "ellipsoidal multivariate normal", 
                 warning("invalid model"))
  return(list(model = model, type = type))
}
getIC <- function(object,criteria){
  if (!is.null(object$models)){
  if(missing(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE) 
  lc <- length(criteria)
  lm <- length(object$models)
  res <- matrix(NA,lm,lc,dimnames=list(1:length(object$models),criteria))
  df <- data.frame(G=1:lm)
  for(i in 1:lm){
    obj <- object$models[[i]]
    if (!is.null(obj$IC)) res[i,] <- sapply(1:lc, FUN=function(j) obj$IC[[criteria[j]]])
    df$G[i] <- obj$G
    if(!is.null(obj$model)) df$model[i] <- as.character(obj$model)
  }
  attributes(res) <- c(attributes(res),df)
  class(res) <- "ContaminatedMixt.IC"
  res
  }
}
print.ContaminatedMixt.IC <- function(x, digits = max(3L, getOption("digits") - 2L),...){
  class(x) <-"matrix"
  d <- data.frame(x)
  d <- data.frame(model= attributes(x)$model,G = attributes(x)$G,d)
  print(d,digits=digits) 
}
getBestModel <- function(object,criterion="BIC",G=NULL,model=NULL){
  #if(!class(object)=="ContaminatedMixt") stop("object is not a class ContaminatedMixt object")
  criterion <- match.arg(criterion,.ICnames())
  n <- whichBest(object=object,criteria=criterion,G=G,model=model)
  if (length(n)>1) stop("More than one model matches the conditions specified.")
  foo <- object$models[[n]]
  object$models <- NULL
  object$models[[1]] <- foo
  invisible(object)
}

whichBest <- function(object, criteria=NULL, G=NULL, model=NULL){
  #returns the best model(s) according to one or more information criteria.
  if(is.null(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE)
  a <- getIC(object,criteria)
  w <- TRUE
  if(!is.null(G))  w <- w & attr(a,"G") %in% G
  if(!is.null(model))  w <- w & attr(a,"model") %in% model
  a <-a[w,,drop=FALSE]
  if (!any(w)) stop("No model matches the conditions specified.")
  best <- strtoi(rownames(a)[sapply(1:length(criteria), function (i) which(t(a[,i])==max(a[,i], na.rm =TRUE)))])
  names(best) <- criteria
  best
}
getCluster <- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$group
}
getPar <- function(object, ...){
  best <- getBestModel(object,...)
  list(prior=best$models[[1]]$prior,mu=best$models[[1]]$mu, Sigma=best$models[[1]]$Sigma,alpha=best$models[[1]]$alpha,eta=best$models[[1]]$eta)
}
getDetection <- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$detection
}
getPosterior<- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$posterior
}
getSize<- function(object, ...){
  best <- getBestModel(object,...)
  table(best$models[[1]]$group)
}