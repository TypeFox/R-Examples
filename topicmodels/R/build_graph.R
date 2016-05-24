## This file is modified from lasso-graph.r which comes with the C
## sources for fitting CTM.

## Graphical model selection using the Lasso, as
## proposed by Meinshausen and Buhlmann

## April, 2007 -- Dave Blei and John Lafferty
##
## To apply this to topic graphs, we take the variational means
## (lambda) for each document, and treat these as data.  We then
## regress each variable (topic) onto the others using the lasso, and
## consider the indices of the non-zero entries as estimates of the
## neighbors of the node in the inverse covariance.  The graph is then
## formed by including an edge if either/both (OR/AND) of the endpoints
## include it in the corresponding penalized regression.

## it's possible to use the lars package as well, with some minor mods

## Inputs
##   file:   n x p data matrix -- e.g., the variational means ("final-lambda.dat")
##   lambda: relative bound on the l1-norm of the parameters, in [0,1]
##   and=T:  if and=T/F then the graph is computed by taking the intersction/union of the nbhds
##
## Output
##   Ihat:   matrix of 0/1, with 1 indicating an edge in the graph

build_graph <- function(x, lambda, and = TRUE) {
  if (!is(x, "CTM")) stop("x needs to be a fitted CTM object")
  gamma <- posterior(x)$topics
  x <- log(sweep(gamma, 1, gamma[,ncol(gamma)], "/"))
  x <- (x - rowMeans(x))/apply(x, 1, sd)
  p <- ncol(x)
  Ihat <- Shat <- matrix(FALSE, p, p)
  for (j in 1:p) {
    DATA <- data.frame(y = x[,j], x[,-j])
    out <- lasso2::l1ce(y ~ ., data = DATA,
                sweep.out = ~1, bound = lambda)
    indices <- (1:p)[-j]
    beta <- coef(out)[-1] # skipping the intercept
    nonzero <- indices[beta > 0]
    Shat[j, nonzero] <- TRUE
  }
  diag(Shat) <- TRUE  
  Ihat <- if (and) Shat & t(Shat) else Shat | t(Shat)
  return(Ihat)
}

