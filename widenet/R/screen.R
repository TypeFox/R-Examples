## functions for screening

## screen.cor and screen.ttests are adapted from the
## SuperLearner package written by Eric Polley



## function to use glmnet to screen the variables
## (running it only once, without cross-validation)

## returns an (integer) index vector for the columns of x
## which are chosen. Length should always be equal to num.vars
## unless num.vars is > ncol(x)

## screen.glmnet <- function(x, y, family = c("gaussian", "binomial"),
##                           num.vars = 100,
##                           nlambda = 100,
##                           alpha = 1,
##                           ...) {

##   family = match.arg(family)
  
##   if(num.vars >= ncol(x)) {
##     warning("no need to screen since x has <= num.vars columns\n",
##             "using all columns")
##     return(1:ncol(x))
##   }

##   glmnet.fit <- glmnet(x, y, family,
##                        nlambda = nlambda,
##                        alpha = alpha,
##                        ...)

##   path.num.vars <- apply(glmnet.fit$beta, 2, function(x) sum(x != 0))

##   if(all(path.num.vars < num.vars)) { ## not enough variables selected

##     ## different warnings depending on alpha
    
##     if(nrow(x) < ncol(x) && alpha == 1) {
##       warning("screening mechanism is selecting some variables at random",
##               "\nsince the glmnet model did not select enough variables.",
##               "\nYou may want to try setting an elastic net penalty",
##               "\nfor the screening, i.e. 0 < screen.alpha < 1")
##     } else {
##       warning("screening mechanism is selecting some variables at random since",
##               "\nthe glmnet model did not select enough variables")
##     }

##     ## max number selected by glmnet

##     max.selected <- max(path.num.vars)
## #oprint(max.selected)
    
##     which.max.selected <- which.max(path.num.vars)

##     logical.index <- glmnet.fit$beta[, which.max.selected] != 0

##     ## get random subset of the rest

##     more.random.indices <- sample((1:ncol(x))[!logical.index],
##                                   size = num.vars - max.selected,
##                                   replace = FALSE)

##     ## return ordered integer vector
    
##     return(as.vector(sort(c(which(logical.index), more.random.indices))))
##   }

##   match.vec <- path.num.vars == num.vars
  
##   if(any(match.vec)) { ## there's at least one match
##                        ## return set given by first match

##     ## find first match

##     first.match <- which(match.vec)[1]
    
##     return(as.vector(which(glmnet.fit$beta[, first.match] != 0)))
##   }

##   ## final possibility -- there is a set given in beta with a number of vars
##   ## greater than num.vars which will have to be subset according to its
##   ## coefficients

##   ## find the smallest number of vars greater than num.vars
  
##   min.greater.val <- min(path.num.vars[path.num.vars > num.vars])

##   first.match <- which(path.num.vars == min.greater.val)[1]

##   coeffs <- glmnet.fit$beta[, first.match]

##   ranks <- rank(abs(coeffs), ties.method = "random")

##   return(as.vector(which(ranks > (length(ranks) - num.vars))))
    
## }


## function to screen using ttests (requires genefilter package)
## requires family = "binomial"

screen.ttest <- function(x, y, family = "binomial",
                         #obsWeights, id,
                         num.vars = 100) {

  if(is.matrix(y))
    stop("The ttest screening function cannot handle a matrix y")

  if(family != "binomial")
    stop("family must be binomial for ttest screening")

  if(num.vars >= ncol(x)) {
    warning("no need to screen since x has <= num.vars columns\n",
            "using all columns")
    return(1:ncol(x))
  }

  ## try to load genefilter package
  
  tryCatch(library(genefilter), error = function(e) {

    stop("screening method is ttest, but unable to load package genefilter.\n",
         "Error message was:\n",
         e$message)
  })
  
  listP <- colttests(x = x, fac = as.factor(y), tstatOnly = FALSE)$p.value
  
  return(as.vector(which(rank(listP, ties.method = "random") <= num.vars)))

}


## function to screen using pearson's correlation

screen.cor <- function(x, y, family,
                       #method = 'pearson', ## always use pearson
                       num.vars = 100) {

  if(num.vars >= ncol(x)) {
    warning("no need to screen since x has <= num.vars columns\n",
            "using all columns")
    return(1:ncol(x))
  }

  listP <- apply(x, 2, function(x.col) { 
                         ifelse(var(x.col) <= 0, 1,
                                cor.test(x.col, y = y,
                                         method = "pearson")$p.value)
                       })
  return(as.vector(which(rank(listP, ties.method = "random") <= num.vars)))
}
