FLXdist <- function(formula, k = NULL, model=FLXMRglm(), components, concomitant=FLXPconstant()) 
{
  mycall <- match.call()
  if(is(model, "FLXM")) model <- list(model)
   
  if (length(k)==1) prior <- rep(1/k, k)
  else {
    prior <- k/sum(k)
  }

  concomitant@x <- matrix(c(1, rep(0, ncol(concomitant@coef))[-1]), nrow = 1)
  prior <- as.vector(evalPrior(prior, concomitant))
  
  lf <- length(formula)
  formula1 <- formula
  if(length(formula[[lf]])>1 && deparse(formula[[lf]][[1]]) == "|")
    formula1[[lf]] <- formula[[lf]][[2]]

  for(n in seq(along.with=model)) {
    if(is.null(model[[n]]@formula))
      model[[n]]@formula <- formula1
    else if(length(model[[n]]@formula) == 3 && model[[n]]@formula[[2]] == ".")
      model[[n]]@formula <- model[[n]]@formula[-2]
    model[[n]]@fullformula <- update.formula(formula1, model[[n]]@formula)
  }
  if (missing(components)) stop("no parameter values specified")
  if (length(components) != length(prior)) stop("components not specified correctly")

  comp <- list()
  for (k in seq(along.with=prior)) {
    comp[[k]] <- list()
    if (length(components[[k]]) != length(model))
      stop("components not specified correctly")
    for (n in seq(along.with=model)) {
      comp[[k]][[n]] <- FLXcomponent(model[[n]],
                                     components[[k]][[n]])
    }
  }
  new("FLXdist", formula=formula, call=mycall, concomitant=concomitant,
      prior=prior, k=length(prior), model=model, components=comp)
}

###**********************************************************

setGeneric("FLXcomponent", function(object, ...) standardGeneric("FLXcomponent"))

setMethod("FLXcomponent", signature(object="FLXM"), function(object, components, ...) {
  df <- numeric()
  with(components, eval(object@defineComponent))
})

##<fixme>##
setMethod("FLXcomponent", signature(object="FLXMRglm"), function(object, components, ...) {
  df <- numeric()
  offset <- NULL
  family <- object@family
  with(components, eval(object@defineComponent))
})

###**********************************************************

setMethod("show", "FLXdist",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nPriors:\n")
    names(object@prior) <- paste("Comp.", seq_along(object@prior), sep="")
    print(object@prior)
    cat("\n")
})

###**********************************************************

evalPrior <- function(prior, concomitant) prior

setGeneric("evalPrior", function(prior, concomitant) standardGeneric("evalPrior"))

setMethod("evalPrior", signature(concomitant="FLXPmultinom"), function(prior, concomitant) {
  exps <- exp(concomitant@x %*% concomitant@coef)
  exps/rowSums(exps)
})
