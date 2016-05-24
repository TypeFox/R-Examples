### These functions are for adaptive MCMC, and only for E[b],
### E[Phi | Phi^{obs}], and E[Phi].

### Initial global storages.
my.set.adaptive <- function(nSave,
    n.aa = NULL, b.DrawScale = NULL,
    n.p = NULL, p.DrawScale = NULL,
    n.G = NULL, phi.DrawScale = NULL,
    n.G.pred = NULL, phi.pred.DrawScale = NULL,
    renew.iter = .CF.AC$renew.iter, adaptive = .CF.CT$adaptive[1]){
  ### Check.
  if(adaptive == "none"){
    total.renew <- 1
  } else{
    total.renew <- ceiling(nSave / renew.iter)
  }

  ### Initials for updating acceptance rate.
  .cubfitsEnv$curr.renew <- 1
  .cubfitsEnv$adaptive <- list(b = list(),
                               p = list(),
                               phi = list(), phi.pred = list())
  .cubfitsEnv$DrawScale <- list(b = list(),
                                p = list(),
                                phi = list(), phi.pred = list())

  ### Assign default to .cubfitsEnv as global variables.
  .cubfitsEnv$all.DrawScale <- list(b = NULL, b.prev = NULL,
                                    p = NULL, p.prev = NULL,
                                    phi = NULL, phi.prev = NULL,
                                    phi.pred = NULL, phi.pred.prev = NULL)

  ### For adaptive rate in parameters.
  if(!is.null(n.aa)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$b[[i.renew]] <- rep(0L, n.aa)
      if(length(b.DrawScale) == 1){
        .cubfitsEnv$DrawScale$b[[i.renew]] <- rep(b.DrawScale, n.aa)
      } else if(length(b.DrawScale) == n.aa){
        .cubfitsEnv$DrawScale$b[[i.renew]] <- b.DrawScale
      } else{
        stop("length of b.DrawScale is incorrect.")
      }
    }

    ### Update scaling factors.
    .cubfitsEnv$all.DrawScale$b <- .cubfitsEnv$DrawScale$b[[1]]
    .cubfitsEnv$all.DrawScale$b.prev <- .cubfitsEnv$all.DrawScale$b
  }

  ### For adaptive rate in prior
  if(!is.null(n.p)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$p[[i.renew]] <- rep(0L, n.p)
      if(length(p.DrawScale) == 1){
        .cubfitsEnv$DrawScale$p[[i.renew]] <- rep(p.DrawScale, n.p)
      } else if(length(p.DrawScale) == n.p){
        .cubfitsEnv$DrawScale$p[[i.renew]] <- p.DrawScale
      } else{
        stop("length of p.DrawScale is incorrect.")
      }
    }

    ### Update scaling factors.
    .cubfitsEnv$all.DrawScale$p <- .cubfitsEnv$DrawScale$p[[1]]
    .cubfitsEnv$all.DrawScale$p.prev <- .cubfitsEnv$all.DrawScale$p
  }

  ### For adaptive rate in expectations.
  if(!is.null(n.G)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi[[i.renew]] <- rep(0L, n.G)
      if(length(phi.DrawScale) == 1){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- rep(phi.DrawScale, n.G)
      } else if(length(phi.DrawScale) == n.G){
        .cubfitsEnv$DrawScale$phi[[i.renew]] <- phi.DrawScale
      } else{
        stop("length of phi.DrawScale is incorrect.")
      }
    }

    ### Update scaling factors.
    .cubfitsEnv$all.DrawScale$phi <- .cubfitsEnv$DrawScale$phi[[1]]
    .cubfitsEnv$all.DrawScale$phi.prev <- .cubfitsEnv$all.DrawScale$phi
  }

  ### For adaptive rate in predictions.
  if(!is.null(n.G.pred)){
    for(i.renew in 1:total.renew){
      .cubfitsEnv$adaptive$phi.pred[[i.renew]] <- rep(0L, n.G.pred)

      if(length(phi.pred.DrawScale) == 1){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- rep(phi.pred.DrawScale,
                                                         n.G.pred)
      } else if(length(phi.pred.DrawScale) == n.G.pred){
        .cubfitsEnv$DrawScale$phi.pred[[i.renew]] <- phi.pred.DrawScale
      } else{
        stop("length of phi.pred.DrawScale is incorrect.")
      }
    }

    ### Update scaling factors.
    .cubfitsEnv$all.DrawScale$phi.pred <- .cubfitsEnv$DrawScale$phi.pred[[1]]
    .cubfitsEnv$all.DrawScale$phi.pred.prev <- .cubfitsEnv$all.DrawScale$phi.pred
  }

  invisible()
} # End of my.set.adaptive().


### Updating function based on variable name and current iteration.
my.update.adaptive <- function(var.name, accept){
  .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] <-
    .cubfitsEnv$adaptive[[var.name]][[.cubfitsEnv$curr.renew]] + accept

  invisible()
} # End of my.update.adaptive().

### copy adaptive information.
my.copy.adaptive <- function(){
  .cubfitsEnv$all.DrawScale$b.prev <- .cubfitsEnv$all.DrawScale$b
  .cubfitsEnv$all.DrawScale$p.prev <- .cubfitsEnv$all.DrawScale$p
  .cubfitsEnv$all.DrawScale$phi.prev <- .cubfitsEnv$all.DrawScale$phi
  .cubfitsEnv$all.DrawScale$phi.pred.prev <- .cubfitsEnv$all.DrawScale$phi.pred

  invisible()
} # End of my.copy.adaptive().

