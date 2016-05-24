setMethod("worth", "raschmix", function(object, difficulty = TRUE,
                                        component = NULL){
  ## get parameters
  coef <- parameters(object, which = "item", difficulty = difficulty,
                     component = component)

  ## apply transformation
  worth <- apply(coef, 2, function(x){
      x[is.na(x)] <- 0 ## if aliased item is inculded as NA
      x - mean(x[is.finite(x)])
  })

  ## include non-identified items (if any)
  rval <- matrix(NA, ncol = ncol(worth), nrow = length(object@identified.items))
  rval[object@identified.items == "0/1", ] <- worth
  rval[object@identified.items == "0", ] <- if(!difficulty) -Inf else Inf
  rval[object@identified.items == "1", ] <- if(!difficulty) Inf else -Inf
  rownames(rval) <- names(object@identified.items)
  colnames(rval) <- colnames(worth)

  return(rval)
})

setMethod("itempar", "raschmix", function(object, ref = NULL, alias = TRUE, ...){

  ## extract parameters
  cf <- parameters(object, which = "item", difficulty = TRUE)
  cf[1,] <- 0
  m <- nrow(cf)
  lbs <- gsub("item.", "", rownames(cf), fixed = TRUE)

  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.vector(ref) && is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.vector(ref) && is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else if (is.matrix(ref) && is.numeric(ref)) {
    stopifnot(nrow(ref) == m && ncol(ref) == m)
  } else stop("Argument 'ref' is misspecified (see ?itempar for possible values).")

  ## if not given, specify contrast matrix
  if (is.matrix(ref)) {
    D <- ref
  } else {
    D <- diag(m)
    D[, ref] <- D[, ref] - 1/length(ref)
  }

  ## apply ref
  cf <- apply(cf, 2, function(x) as.vector(D %*% x))
  
  ## items solved by no or all subjects
  rval <- matrix(NA, ncol = ncol(cf), nrow = length(object@identified.items))
  rval[object@identified.items == "0/1", ] <- cf
  rval[object@identified.items == "0", ] <-  Inf
  rval[object@identified.items == "1", ] <- -Inf
  
  rownames(rval) <- names(object@identified.items)
  colnames(rval) <- colnames(cf)

  ## ## if vcov requested: adjust existing vcov
  ## if (vcov) {
  ##   warning("Variance covariance matrix not implemented for Rasch mixture models.")
  ## }
  ## vc <- matrix(NA, nrow = length(object@identified.items), ncol = length(object@identified.items))
  
  ## ## set labels
  ## rownames(rval) <- rownames(vc) <- colnames(vc) <- names(object@identified.items) #lbs


  ## process argument alias
  if (!alias) {
      if (is.matrix(ref)) {
          ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
          stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
      } else {
          aliased <- which(object@identified.items == "0/1")[ref[1]]
          rval <- rval[-aliased,]
          #vc <- vc[-aliased, -aliased]
          alias <- paste0("I", aliased)
          names(alias) <- names(object@identified.items)[aliased]
      }
  }
  
  ## ## setup and return result object
  ## rv <- structure(rval, class = "itempar", model = "raschmix", ref = ref, alias = alias, vcov = vc)
  ## return(rv)
  return(rval)
    
})


scoreProbs <- function(object, component = NULL, simplify = TRUE){

  ## get score parameters
  delta <- parameters(object, which = "score", component = component,
                      simplify = FALSE)
  delta <- lapply(delta, function(x) x$score)

  ## set up design matrix
  m <- sum(object@identified.items == "0/1")
  rs <- which(object@rawScoresData > 0)
  nscores <- length(rs)
  switch(object@scores,
    "saturated" = {
      xaux <- diag(nscores)
      extra <- 0
    },
    "meanvar" = {
      rr <- 1:(m-1L)
      xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
      extra <- NULL
    },      
    "constant" = {
      xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
      extra <- 0
    }
    )

  ## calculate score probabilities
  scores <- sapply(delta, simplify = simplify, FUN = function(delta){

    ## in case any of the scores was not estimated in this component
    ## (due to weight 0)
    if(object@scores == "saturated"){
        ref <- which(is.na(delta))
        delta <- delta[is.finite(delta)]
        xaux <- diag(length(delta) + 1)
    }    

    ## probabilities (colSums instead of %*% to get na.rm = TRUE)
    eta <- colSums(c(extra, delta) * t(xaux), na.rm = TRUE)
    psi <- exp(eta) / sum(exp(eta))
    
    ## check for remaining problems
    if(any(is.na(psi))) stop("some parameters in score model not identified")

    ## for saturated model:
    ## return probability 0 for scores not present in the data
    if (object@scores == "saturated"){
      psi.full <- numeric(m - 1L)
      psi.full[c(ref, as.numeric(names(delta)))] <- psi
      psi <- psi.full
    }

    # conditional MLE -> r = 0 and r = "number of items" are excluded 
    psi <- psi * (1 - sum(object@extremeScoreProbs))
    psi <- c(object@extremeScoreProbs[1], psi, object@extremeScoreProbs[2])
    
    return(psi)
  })

  return(scores)
}



setMethod("parameters", "raschmix", function(object,
       which = c("model", "item", "score", "concomitant"),
       difficulty = TRUE, component = NULL, simplify = TRUE){
  
  which <- match.arg(which)
  which.flx <- if (which %in% c("item", "score")) "model" else which

  if (is.null(component)) component <- 1:object@k
  drop <- length(component) > 1
  
  ## call flexmix method
  para <- callNextMethod(object, component = component, model = NULL,
                         which = which.flx, simplify = FALSE, drop = drop)

  if (which != "concomitant"){

    if (length(component) == 1) para <- para[[1]]

    para <- lapply(para, function(comp){

#      comp$item <- comp$item[is.finite(comp$item)]
#      comp$item <- comp$item[-1]
      comp$item[min(which(is.finite(comp$item)))] <- NA  
      if (!difficulty) comp$item <- -comp$item

#      alias <- if(length(comp$score) > 0L) !is.finite(comp$score) else logical(0L)
#      comp$score <- comp$score[!alias]
      if(object@scores == "saturated"){
        score.full <- rep(NA, length.out = (length(comp$score)+1))
        names(score.full) <- seq_along(score.full)
        score.full[as.numeric(names(comp$score))] <- comp$score
        comp$score <- score.full
      }
      if (object@scores == "meanvar"){
        names(comp$score) <- c("location", "dispersion")
      }
      
      if (which == "score") comp$item <- NULL
      if (which == "item") comp$score <- NULL
      return(comp)
    })

    if (simplify){
      para <- lapply(para, unlist)
      para <- do.call("cbind", para)
    }
    ## <FIXME> for object@scores == "constant" & which == "score":
    ## return sth "empty" rather than NULL? </FIXME>

  }
  return(para)
})



## <INFO>
## summary-method for raschmix objects:
## it shows:
## call (to raschmix), information about the mixture, item parameters, logLik, AIC, and BIC
## it doesn't show (yet):
## standard errors for item parameters (needs refit function), coefficients for concomitant model
## other notes:
## flexmix separates display of mixture from display of components (summary vs parameters methods) --> stick with that?
## an option to select a certain component is not implemented as it's a summary of the mixture
## there is no summary method for stepFlexmix objects, just a show method -- this gets updated to use
## the raschmix call instead of the flexmix call but is otherwise left untouched
## </INFO>

## modified code from flexmix
setClass("summary.raschmix",
         representation(itemParaTab = "ANY"),
         contains = "summary.flexmix")

setMethod("show", "summary.raschmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    print(object@comptab, digits=3)
    cat("\nItem Parameters:\n")
    print(object@itemParaTab)
    cat("\n")
    print(object@logLik)
    cat("AIC:", object@AIC, "  BIC:", object@BIC, "\n")
    cat("\n")    
})

setMethod("summary", "raschmix",
function(object, eps=1e-4, ...){    
    z <- new("summary.raschmix",
             call = object@call,
             AIC = AIC(object),
             BIC = BIC(object),
             logLik = logLik(object))

    TAB <- data.frame(prior=object@prior,
                      size=object@size)
    rownames(TAB) <- paste("Comp.", seq_len(nrow(TAB)), sep="")
    TAB[["post>0"]] <- colSums(object@posterior$scaled > eps)
    TAB[["ratio"]] <- TAB[["size"]]/TAB[["post>0"]]
    
    z@comptab = TAB
    z@itemParaTab <- worth(object)
    z
})


## <INFO>
## logLik-method for "stepRaschmix" doesn't display df correctly,
## but logLik in not indented for this use anyway
##
## refit
## certainly for raschmix, also for "stepRaschmix"?
##
## ## weights
## ## for raschmix
## setGeneric("weights")
## setMethod("weights", "raschmix", function(object){
##   if(!is.null(object@weights)) object@weights else rep(1, object@nobs)
## })
##
## ## EIC for "stepRaschmix"?
## </INFO>
