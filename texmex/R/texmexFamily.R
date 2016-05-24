texmexFamily <-
    # Create an object of class 'texmexFamily'. It is allowable to have
    # info, start and resid as NULL, but all other information must
    # be provided by the user.
function(name, log.lik, param, info=NULL, start=NULL, resid=NULL,
                         rl, delta, endpoint, density, rng, prob, quant){
    res <- list(name=name, log.lik=log.lik, param=param, info=info, start=start,
                resid=resid, rl=rl, delta=delta, endpoint=endpoint, density=density,
                rng=rng, prob=prob, quant=quant)

    oldClass(res) <- 'texmexFamily'
    res
}


print.texmexFamily <- function(x, ...){
    cat('Family:      ', x$name, '\n')

    invisible()
}

summary.texmexFamily <- function(object, ...){
  if (is.null(object$info)){ info <- 'Numerical approximation' }
  else { info <- 'Closed form' }

  print.texmexFamily(object, verbose=FALSE, ...)
  cat('Parameters:  ', object$param, '\n')
  cat('Information: ', info, '\n')
  invisible()
}