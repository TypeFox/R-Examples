# Plots two dimensional graph of likelihood
plotLikelihood.2d <- function(hypothesis, which=c(1, 2), large=100, N=20,
                              arguments=NULL, x=NULL, y=NULL,
                              logObjective=TRUE, logDegradation=TRUE,
                              contours=list(), ...) {
  # define ggplot variables to avoid NOTE from CRAN
aes <- xlab <- ylab <- geom_tile <- labs <- stat_contour <- NULL

  if(length(which) != 2)
    stop("Argument 'which' to plotLikelihood.2d' should be of length 2")

  sanity.check = getFromNamespace("sanity.check", "likeLTD")
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  # Create objective function
  if(logObjective) creator <- create.likelihood.log

  params = optimisation.params( hypothesis, verbose=FALSE,
                                logObjective=logObjective,
                                logDegradation=logDegradation, 
                                arguments=arguments, zero=0e0, ... )
  if(any(which > length(params$par))) stop("'which' argument it too large")
  if(any(which < 1)) stop("'which' argument it too small")

  # Now figure out extents.
  xRange = c(params$lower[which[[1]]], params$upper[which[[1]]]) 
  yRange = c(params$lower[which[[2]]], params$upper[which[[2]]]) 

  # Replace infinities with large number.
  if(length(large) == 1) large = rep(large, 2)
  if(abs(xRange[[1]]) == Inf) xRange[[1]] = sign(xRange[[1]]) * large[[1]]
  if(abs(xRange[[2]]) == Inf) xRange[[2]] = sign(xRange[[2]]) * large[[1]]
  if(abs(yRange[[1]]) == Inf) yRange[[1]] = sign(yRange[[1]]) * large[[2]]
  if(abs(yRange[[2]]) == Inf) yRange[[2]] = sign(yRange[[2]]) * large[[2]]
  
  # Now create map
  if(length(N) == 1) N = rep(N, 2)
  if(is.null(x))
    x = ((1:N[[1]]) - 1) / (N[[1]]-1) * (xRange[[2]] - xRange[[1]]) + xRange[[1]]
  if(is.null(y)) 
    y = ((1:N[[2]]) - 1) / (N[[2]]-1) * (yRange[[2]] - yRange[[1]]) + yRange[[1]]
  

  map = matrix(0, nrow=length(x), ncol=length(y))
  for(i in 1:length(x)) for(j in 1:length(y)) {
    newargs = params$par
    newargs[[which[[1]]]] = x[[i]]
    newargs[[which[[2]]]] = y[[j]]
    map[i, j] = params$fn(newargs)
  }

  amap = data.frame(x=rep(x, length(y)), y=rep(y, rep(length(x), length(y))),
                    z=c(map))
  z = amap$z # Avoids buggy warning when checking package.
  title = sprintf("nUnknowns=%d, dropin=%d", hypothesis$nUnknowns,
                  hypothesis$doDropin)
  fill = "Likelihood"
  if(logObjective) fill = "Log-likelihood"
  ggplot(amap, aes(x=x, y=y, z=z))                  +
     xlab(names(unlist(arguments))[[which[[1]]]])   +
     ylab(names(unlist(arguments))[[which[[2]]]])   +
     geom_tile(aes(fill=z))                         +
     labs(fill=fill)                                + 
     do.call(stat_contour, contours)
}
