## Find the indices into inla output config structures corresponding
## to a specific predictor, effect name, or inla.stack tag.
##
## result : an inla object
inla.output.indices = function(result, name=NULL, stack=NULL, tag=NULL, ...)
{
    if (!is.null(name) && !is.null(tag)) {
        stop("At most one of 'name' and 'tag' may be non-null.")
    }
    if (!is.null(tag)) {
        if (is.null(stack)) {
            stop("'tag' specified but 'stack' is NULL.")
        }
        tags <- names(stack$data$index)
        if (!(tag %in% tags)) {
            stop("'tag' not found in 'stack'.")
        }
    } else if (is.null(name) || (!is.null(name) && (name == ""))) {
        if ("APredictor" %in% result$misc$configs$contents$tag) {
            name <- "APredictor"
        } else {
            name <- "Predictor"
        }
    }

    ## Find variables
    if (!is.null(name)) {
        if (!(name %in% result$misc$configs$contents$tag)) {
          stop("'name' not found in result.")
        }
        nameindex <- which(result$misc$configs$contents$tag == name)
        index <- (result$misc$configs$contents$start[nameindex] - 1L +
                  seq_len(result$misc$configs$contents$length[nameindex]))
    } else { ## Have tag
        index <- stack$data$index[[tag]]
    }

    index
}

private.simconf.link <- function(res,links,trans=TRUE)
{
  if(trans){
    n = length(res$a)
    res$a.marginal = sapply(1:n, function(i) private.link.function(
                                      res$a.marginal[i],links[i],inv=TRUE))
    res$b.marginal = sapply(1:n, function(i) private.link.function(
                                      res$b.marginal[i],links[i],inv=TRUE))
    res$a = sapply(1:n, function(i) private.link.function(
                                      res$a[i],links[i],inv=TRUE))
    res$b = sapply(1:n, function(i) private.link.function(
                                      res$b[i],links[i],inv=TRUE))
  }
  return(res)
}

private.link.function <- function(x, link, inv=FALSE)
{
  if (is.na(link)) {
    link = "identity"
  }
  return(do.call(paste("inla.link.", link, sep=""),list(x=x, inv=inv)))
}

private.get.config <- function(result,i)
{
  mu=result$misc$configs$config[[i]]$mean
	Q=forceSymmetric(result$misc$configs$config[[i]]$Q)
	vars = diag(result$misc$configs$config[[i]]$Qinv)
	lp = result$misc$configs$config[[i]]$log.posterior
  return(list(mu=mu,Q=Q,vars=vars,lp=lp))
}

## Calculate the marginal probability for X_i>u or X_i<u.
## Note that the index 'i' refers to a location in the linear
## predictor if predictor==TRUE, whereas it refers to a location
## in the random effect vector otherwise.
inla.get.marginal <- function(i, u,result,effect.name=NULL, u.link, type)
{
  link = result$misc$linkfunctions$names[result$misc$linkfunctions$link][i]

  if(is.null(effect.name) && u.link == TRUE){
    #calculate marginals using fitted values
    #if(link=="identity" || is.na(link)){
      marg.p = result$marginals.fitted.values[[i]]
    #} else {

      #marg.p = INLA::inla.tmarginal(function(x)
      #                        private.link.function(x,link,inv=TRUE),
      #                            result$marginals.fitted.values[[i]])
    #}
  } else if (is.null(effect.name)){
      # Calculate marginals using linear predictor
      marg.p = result$marginals.linear.predictor[[i]]
  } else {
      # Calculate marginals using a random effect
      marg.p = result$marginals.random[[effect.name]][[i]]
  }

	  if(type=='<'){
		  return(INLA::inla.pmarginal(u,marg.p))
	  } else {
		  return(1-INLA::inla.pmarginal(u,marg.p))
	  }
  }