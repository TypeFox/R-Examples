
"prior.glm.control" <- 
  function(beta.prior = c("flat", "normal", "fixed"), beta = NULL, beta.var.std = NULL, 
           sigmasq.prior = c("uniform", "sc.inv.chisq", "reciprocal", "fixed"), sigmasq = NULL, df.sigmasq = NULL,
           phi.prior = c("uniform","exponential", "fixed", "squared.reciprocal", "reciprocal"), phi = NULL, 
           phi.discrete = NULL, tausq.rel = 0)
{
  beta.prior <- match.arg(beta.prior)
  if(beta.prior == "fixed" & is.null(beta))
    stop("argument \"beta\" must be provided with fixed value for this parameter")
  if(beta.prior == "normal"){
    if(is.null(beta) | is.null(beta.var.std))
      stop("arguments \"beta\" and \"beta.var.std\" must be provided when using normal prior for the parameter beta")
    if((length(beta))^2 != length(beta.var.std))
      stop(" beta and beta.var.std have incompatible dimensions")
    if(any(beta.var.std != t(beta.var.std)))
      stop(" non symmetric matrix in beta.var.std")
    if(inherits(try(chol(beta.var.std)), "try-error"))
      stop(" matrix in beta.var.std is not positive definit")
  }
  ##
  sigmasq.prior <- match.arg(sigmasq.prior)
  if(sigmasq.prior == "fixed" & is.null(sigmasq))
    stop("argument \"sigmasq\" must be provided when the parameter sigmaq is fixed")
  if(sigmasq.prior == "sc.inv.chisq")
    if(is.null(sigmasq) | is.null(df.sigmasq))
      stop("arguments \"sigmasq\" and \"df.sigmasq\" must be provided for inverse chisq prior")
  if(!is.null(sigmasq))
    if(sigmasq < 0) stop("negative values not allowed for \"sigmasq\"")
  if(sigmasq.prior == "reciprocal"){
    warning("This choice of sigmasq.prior gives an improper posterior !!!!!!! \n")
    sigmasq <- 0
    df.sigmasq <- 0
  }
  if(sigmasq.prior == "uniform"){
    sigmasq <- 0
    df.sigmasq <- -2
  }
  ##
  if(!is.null(phi) && length(phi) > 1)
    stop("prior.glm.control: length of phi must be one. ")
  if(is.numeric(phi.prior)){
    phi.prior.probs <- phi.prior
    phi.prior <- "user"
    if(is.null(phi.discrete))
      stop("prior.glm.control: argument phi.discrete with support points for phi must be provided\n")
    if(length(phi.prior.probs) != length(phi.discrete))
      stop("prior.glm.control: user provided phi.prior and phi.discrete have incompatible dimensions\n")
    if(round(sum(phi.prior.probs), digits=8) != 1)
      stop("prior.glm.control: prior probabilities provided for phi do not sum up to 1")
  }
  else phi.prior <- match.arg(phi.prior)
  if(phi.prior == "fixed"){
    if(is.null(phi)){
      stop("argument \"phi\" must be provided with fixed prior for this parameter")
    }
    phi.discrete <- phi
  }
  else{
    if(phi.prior == "exponential" & (is.null(phi) | (length(phi) > 1)))
      stop("argument \"phi\" must be provided when using the exponential prior for the parameter phi")
    if(!is.null(phi.discrete)){
      discrete.diff <- diff(phi.discrete)
      if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
        stop("The current implementation requires equally spaced values in the argument \"phi.discrete\"\n")
    }
    if(phi.prior != "exponential") phi <- NULL
    if(sigmasq.prior == "fixed") stop("option for fixed sigmasq and random phi not implemented")
  }
  if(any(phi.discrete < 0))
    stop("negative values not allowed for parameter phi")
  ##
  if(is.null(tausq.rel)) stop("argument \"tausq.rel\" must be provided")
  ##
  ip <- list(beta=list(), sigmasq=list(), phi=list())
  ##
  if(beta.prior == "fixed"){
    ip$beta$status <- "fixed"
    ip$beta$fixed.value <- beta 
  }
  else{
    ip$beta <- list(dist = beta.prior)
    if(beta.prior == "flat")
      ip$beta$pars <- c(0, +Inf)
    if(beta.prior == "normal"){
      if(length(beta) == 1)
        ip$beta$pars <- c(mean=beta, var.std=beta.var.std)
      else
        ip$beta$pars <- list(mean=beta, var.std=beta.var.std)
    }
  }
  ##
  if(sigmasq.prior == "fixed"){
    ip$sigmasq$status <- "fixed"
    ip$sigmasq$fixed.value <- sigmasq 
  }
  else{
    ip$sigmasq <- list(dist = sigmasq.prior)
    if(sigmasq.prior == "reciprocal")
      ip$sigmasq$pars <- c(df=0, var=+Inf)
    if(sigmasq.prior == "uniform")
      ip$sigmasq$pars <- c(df=-2, var=+Inf)
    if(sigmasq.prior == "sc.inv.chisq")
      ip$sigmasq$pars <- c(df=df.sigmasq, var=sigmasq)
  }
  ##
  if(phi.prior == "fixed"){
    ip$phi$status <- "fixed"
    ip$phi$fixed.value <- phi
  }
  else{
    ip$phi$dist <- phi.prior
    if(is.null(phi.discrete))
      stop("phi.discrete must be given when parameter phi is random")
    else{
      pd <- as.vector(phi.discrete)
      names(pd) <- NULL
      ip$phi$probs <- switch(phi.prior,
                             uniform = rep(1, length(pd)),
                             exponential = dexp(pd, rate=1/phi),
                             squared.reciprocal = ifelse((pd > 0), 1/(pd^2),0),
                             reciprocal = ifelse((pd > 0), 1/pd, 0),
                             user = phi.prior.probs)
      names(ip$phi$probs) <- phi.discrete
    }
    if(phi.prior == "exponential") ip$phi$pars <- c(ip$phi$pars, exp.par=phi)
    ip$phi$probs <- ip$phi$probs/sum(ip$phi$probs)
  }
  ##
  ip$tausq.rel <- list(status = "fixed", fixed.value = tausq.rel)
  ##
  res <- list(beta.prior = beta.prior, beta = beta, beta.var.std = beta.var.std, sigmasq.prior = sigmasq.prior, sigmasq = sigmasq, 
              df.sigmasq = df.sigmasq, phi.prior = phi.prior, phi = phi, 
              phi.discrete = phi.discrete, tausq.rel = tausq.rel, priors.info = ip)
  class(res) <- "prior.geoRglm"
  return(res)
}


".prior.glm.check.aux" <-
  function(prior, fct)
{
  if(class(prior) != "prior.geoRglm"){
    if(!is.list(prior))
      stop(paste(fct,": argument prior only takes a list or an output of the function prior.glm.control"))
    else{
      prior.names <- c("beta.prior", "beta", "beta.var.std", "sigmasq.prior",
                       "sigmasq", "df.sigmasq", "phi.prior", "phi", "phi.discrete", "tausq.rel") 
      prior <- .object.match.names(prior,prior.names)
      if(is.null(prior$beta.prior)) prior$beta.prior <- "flat"
      if(is.null(prior$sigmasq.prior)) prior$sigmasq.prior <- "uniform"
      if(is.null(prior$phi.prior)) prior$phi.prior <- "uniform"
      if(is.null(prior$tausq.rel)) prior$tausq.rel <- 0
      prior <- prior.glm.control(beta.prior = prior$beta.prior,
                                 beta = prior$beta, beta.var.std = prior$beta.var.std,
                                 sigmasq.prior = prior$sigmasq.prior,
                                 sigmasq = prior$sigmasq,  df.sigmasq = prior$df.sigmasq,
                                 phi.prior = prior$phi.prior,
                                 phi = prior$phi, phi.discrete = prior$phi.discrete, 
                                 tausq.rel = prior$tausq.rel)
    }
  }
  return(prior)
}


"image.glm.krige.bayes" <-
  function (x, locations, borders, 
            values.to.plot = c("median", "uncertainty",
              "quantiles", "probabilities", "simulation"),
            number.col, coords.data,
            x.leg, y.leg, messages, ...) 
{
  ## apart from the arguments for "values.to.plot", this function is identical to "image.krige.bayes"
  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
##  attach(x, pos=2, warn.conflicts=FALSE)
##  on.exit(detach("x"))
  if(missing(locations))
    locations <-  eval(attr(x, "prediction.locations"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(mode(values.to.plot) != "numeric")
    values.to.plot <- match.arg(values.to.plot, choices = c("median", "uncertainty","quantiles", "probabilities", "simulation"))
  if(missing(borders)){
    if(!is.null(attr(x, "borders"))) borders.arg <- borders <- eval(attr(x, "borders"))
    else borders.arg <- borders <- NULL
  }
  else{
    borders.arg <- borders
    if(is.null(borders)) borders <- eval(attr(x, "borders"))
  }
  if(missing(number.col)) number.col <- NULL
  if(missing(coords.data)) coords.data <- NULL
  else
    if(coords.data) coords.data <-  eval(attr(x, "data.locations"))
  if(missing(x.leg)) x.leg <- NULL
  if(missing(y.leg)) y.leg <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D')
    do.call("plot.1d", c(list(x = values.to.plot,
                              x1vals = unique(round(locations[,1], digits=12))), .ldots.set(ldots, type="plot.1d",
                                   data="prediction")))
  else{
    ldots.image <- .ldots.set(ldots, type="image", data="prediction")
    locations <- .prepare.graph.krige.bayes(obj=x,
                                           locations=locations,
                                           borders=borders,
                                           borders.obj = eval(attr(x, "borders")),
                                           values.to.plot=values.to.plot,
                                           number.col = number.col,
                                           xlim= ldots.image$xlim,
                                           ylim= ldots.image$ylim,
                                           messages=messages, bound=TRUE)
    do.call("image", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.image))
    if(!is.null(coords.data)) points(coords.data)
    if(!is.null(borders.arg)) polygon(borders, lwd=2)
    if(is.null(ldots$col)) ldots$col <- heat.colors(12)
    if(!is.null(x.leg) & !is.null(y.leg)){
      do.call("legend.krige", c(list(x.leg=x.leg, y.leg=y.leg,
                                     values=locations$values),
                                     ldots))
    }
  }
  return(invisible())  
}


"persp.glm.krige.bayes" <-
  function (x, locations, borders, 
            values.to.plot = c("median", "uncertainty",
              "quantiles", "probabilities", "simulation"), number.col, messages, ...) 
{
  ## apart from the "values.to.plot" argument, this function is identical to "image.glm.krige.bayes"
  ldots <- match.call(expand.dots = FALSE)$...
  if(missing(x)) x <- NULL
  ## attach(x, pos=2, warn.conflicts=FALSE)
  ## on.exit(detach("x"))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2) stop("locations must be a matrix or data-frame with two columns")
  if(mode(values.to.plot) != "numeric"){
    values.to.plot <- match.arg(values.to.plot,
                                choices = c("median", "uncertainty",
                                  "quantiles", "probabilities", "simulation"))
  }
  if(missing(borders)) borders <- NULL
  if(missing(number.col)) number.col <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D')
    do.call("plot.1d", c(list(x = values.to.plot,
                              x1vals = unique(round(locations[,1], digits=12))), .ldots.set(ldots, type="plot.1d",
                                    data="prediction")))
  else{
    ldots.persp <- .ldots.set(ldots, type="persp", data="prediction")
    locations <- .prepare.graph.krige.bayes(obj=x, locations=locations,
                                           borders=borders,
                                           borders.obj = eval(attr(x, "borders")),
                                           values.to.plot=values.to.plot,
                                           xlim= ldots.persp$xlim,
                                           ylim= ldots.persp$ylim,
                                           number.col = number.col,
                                           messages=messages, bound=TRUE)
    do.call("persp", c(list(x=locations$x, y=locations$y,
                            z=locations$values), ldots.persp))
  }
  return(invisible())
}




### The cases where some of the parameters are fixed are not implemented yet.

"hist.glm.krige.bayes" <-
  function(x, pars, density.est = TRUE,
           histogram = TRUE, ...)
{
  ## we need to construct the object such that it fits into geoR function hist.krige.bayes
  kb <- list(posterior=list(),prior=list())
  kb$prior$tausq.rel <- list(status="fixed") 
  kb$prior$phi <- list(status=ifelse(!is.null(x$prior$phi$status), x$prior$phi$status, "random"))
  kb$prior$sigmasq <- list(status=ifelse(!is.null(x$prior$sigmasq$status),x$prior$sigmasq$status, "random"))
  if(is.vector(x$posterior$beta$sample)){
    n.sim <- length(x$posterior$beta$sample)
    if(kb$prior$phi$status =="random"){
      kb$posterior$sample <- as.data.frame(cbind(x$posterior$beta$sample,
                                                 x$posterior$sigmasq$sample, x$posterior$phi$sample)) 
    }
    else{
      if(kb$prior$sigmasq$status =="random"){
        kb$posterior$sample <- as.data.frame(cbind(x$posterior$beta$sample,
                                                   x$posterior$sigmasq$sample, rep(9999,n.sim))) 
      }
      else kb$posterior$sample <- as.data.frame(cbind(x$posterior$beta$sample,
                                                      rep(-99,n.sim), rep(9999,n.sim))) 
    }
    names(kb$posterior$sample) <- c("beta", "sigmasq", "phi")
  }
  else{
    n.sim <- nrow(x$posterior$beta$sample)
    if(kb$prior$phi$status =="random"){
      kb$posterior$sample <- as.data.frame(cbind(t(x$posterior$beta$sample),
                                                 x$posterior$sigmasq$sample, x$posterior$phi$sample))
    }
    else{
      if(kb$prior$sigmasq$status =="random"){
        kb$posterior$sample <- as.data.frame(cbind(t(x$posterior$beta$sample),
                                                   x$posterior$sigmasq$sample, rep(9999,n.sim)))
      }
      else kb$posterior$sample <- as.data.frame(cbind(t(x$posterior$beta$sample),
                                                      rep(-99,n.sim), rep(9999,n.sim)))
    }
    names(kb$posterior$sample) <- c(names(x$posterior$beta$mean), "sigmasq", "phi")
  }
  kb$posterior$sample$tausq.rel <- rep(-99,n.sim)
  hist.krige.bayes(x=kb, pars=pars, density.est = density.est, histogram = histogram)
  return(invisible())
}
