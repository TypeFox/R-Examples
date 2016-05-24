getIC <- function(object,criteria){
  if(missing(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE) 
  lc <- length(criteria)
  lm <- length(object$models)
  res <- matrix(NA,lm,lc,dimnames=list(1:length(object$models),criteria))
  df <- data.frame(k=1:lm)
  for(i in 1:lm){
    obj <- object$models[[i]]
    res[i,] <- sapply(1:lc, FUN=function(j) obj$IC[[criteria[j]]])
    df$k[i] <- obj$k
    fam <- .getFamily(object,i)
    if(!is.null(fam)) df$familyY[i] <-.getFamily(object,i)
    if(!is.null(obj$concomitant$normal.model)) df$normal.model[i] <- as.character(obj$concomitant$normal.model)
  }
  attributes(res) <- c(attributes(res),df)
  class(res) <- "cwm.IC"
  res
}
.getFamily <- function (object,i){
  if (length(object$models[[i]]$GLModel)>0){
    fam <- family(object$models[[i]]$GLModel[[1]]$model)   
    paste0(fam$family, "(",fam$link,")")}
  } 
print.cwm.IC <- function(x, digits = max(3L, getOption("digits") - 2L),...){
  class(x) <-"matrix"
  d <- data.frame(x)
  if (!is.null(attributes(x)$normal.model)){
    d <- data.frame(model= attributes(x)$normal.model,d)
  }
  if (!is.null(attributes(x)$familyY)){
    d <- data.frame(familyY= attributes(x)$familyY,d)
  }
  d <- data.frame(k = attributes(x)$k,d)
  
  print(d,digits=digits) 
}
whichBest <- function(object, criteria=NULL, k=NULL, modelXnorm=NULL, familyY=NULL){
#returns the best model(s) according to one or more information criteria.
  if(is.null(criteria)) criteria <- .ICnames()
  else criteria <- match.arg(criteria, .ICnames(),several.ok =TRUE)
  
  if (!is.null(familyY)){
    familyY <- .familyY(familyY)
    familyY <- sapply(1:length(familyY),function(i) paste0(familyY[[i]]$family,"(",familyY[[i]]$link,")"))
  }
  a   <- getIC(object,criteria)
  w <- TRUE
  if(!is.null(k))  w <- w & attr(a,"k") %in% k
  if(!is.null(modelXnorm))  w <- w & attr(a,"normal.model") %in% modelXnorm
  if(!is.null(familyY))  w <- w & attr(a,"familyY") %in% familyY
  a <-a[w,,drop=FALSE]
  if (!any(w)) stop("No model matches the conditions specified.")
  best <- strtoi(rownames(a)[sapply(1:length(criteria), function (i) which(t(a[,i])==max(a[,i])))])
  names(best) <- criteria
  best
}
getBestModel <- function(object,criterion="BIC",k=NULL,modelXnorm=NULL, familyY=NULL){
  if(!class(object)=="cwm") stop("object is not a class cwm object")
  criterion <- match.arg(criterion,.ICnames())
  n <- whichBest(object=object,criteria=criterion,k=k,modelXnorm=modelXnorm,familyY=familyY)
  if (length(n)>1) stop("More than one model matches the conditions specified.")
  foo <- object$models[[n]]
  object$models <- NULL
  object$models[[1]] <- foo
  invisible(object)
}
getCluster <- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$cluster
}
getPosterior<- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$posterior
}
getSize<- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$size
}
getParPrior<- function(object, ...){
  best <- getBestModel(object,...)
  best$models[[1]]$prior
}
getParConcomitant<- function(object,name=NULL, ...){
  m <- getBestModel(object,...)$models[[1]]$concomitant
  if (!is.null(name)) name <- match.arg(name,c("multinomial","poisson","normal","binomial"),several.ok = TRUE)
  else name <- c("multinomial","poisson","normal","binomial") 
  res <- list()
  if ("normal" %in% name & !is.null(m[["normal.mu"]])) res <- c(res,m["normal.mu"], m["normal.Sigma"])
  if ("binomial" %in% name & !is.null(m[["binomial.p"]])) res <- c(res, m["binomial.p"])
  if ("poisson" %in% name & !is.null(m[["poisson.lambda"]])) res <- c(res, m["poisson.lambda"])
  if ("multinomial" %in% name & !is.null(m[["multinomial.probs"]])) res <- c(res, m["multinomial.probs"])
  res
}
getSize<- function(object, ...){
  getBestModel(object,...)$models[[1]]$size
}
getParXmult <- function(object, ...)getParConcomitant(object, name="multinomial",...)
getParXpois <- function(object, ...)getParConcomitant(object, name="poisson",...)
getParXnorm <- function(object, ...)getParConcomitant(object, name="normal",...)
getParXbin <- function(object, ...)getParConcomitant(object, name="binomial",...)
.mygrep <- function(string,name) {
  substring(string, 1, 2) %in% substring(name,1,2)
}
getPar <- function(object, ...){
  best <- getBestModel(object,...)
  c(prior=getParPrior(best),getParGLM(best),concomitant=getParConcomitant(best))
}
getParGLM <- function(object, ...){
  best <- getBestModel(object,...)
  obj  <- best$models[[1]]
  if (!is.null(obj$GLModel)){
  lr <- lapply(seq_len(obj$k), function(i){
    par <- obj$GLModel[[i]]
    c(list(coefficients=par$model$coefficients),par[-1])
  })
  names(lr) <- paste0("GLMComp.",seq_len(obj$k))
  lr
  } else NULL
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
.ICnames <- function(x=NULL){
  v <- c("AIC", "AICc", "AICu", "AIC3", "AWE","BIC", "CAIC", "ICL")
  if (is.null(x)) return(v)
  else return (v[x])
}


.familyY <- function(familyY){
  if (!is.list(familyY) || !is.null(familyY$family)) familyY <- list(familyY=familyY)
  
  for (i in 1:length(familyY))
    if (is.character(familyY[[i]])){
      familyYname <- match.arg(familyY[[i]],c("gaussian","poisson","binomial","Gamma","student.t","inverse.gaussian"))
      familyY[[i]] <- do.call(familyYname, list())
    } 
  else if (is.function(familyY[[i]])) familyY[[i]] <- familyY[[i]]()
  familyY
}