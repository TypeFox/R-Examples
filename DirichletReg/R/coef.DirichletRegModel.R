coef.DirichletRegModel <- function(object, type = c("both", "beta", "gamma"), ...){

  type <- match.arg(type, c("both", "beta", "gamma"))

  cc <- object$coefficients
  names(cc) <- object$coefnames

  if(object$parametrization == "common"){

    ind <- cumsum(c(1L, object$n.vars))
    cc <- lapply(seq_along(ind)[-1L], function(i) cc[seq.int(from = ind[i - 1L], to = ind[i] - 1L)])
    names(cc) <- object$varnames
    return(cc)

  } else {

    ind <- cumsum(c(1L, object$n.vars))
    bb <- list()

    bb_par <- cc[-seq.int(from = length(cc) - ncol(object$Z) + 1L, to = length(cc))]
    iter <- 0L
    for(i in seq_len(object$dims)){
      if(i == object$base){
        bb[i] <- list(NULL)
      } else {
        iter <- iter + 1L
        bb[[i]] <- bb_par[seq.int(from = ind[iter], to = ind[iter + 1L] - 1L)]
      }
    }
    names(bb) <- c(object$varnames)
    gg <- list("gamma" = cc[seq.int(from = length(cc) - ncol(object$Z) + 1L, to = length(cc))])

    if(type == "both"){
      return(list("beta"=bb, "gamma"=gg))
    } else if(type == "beta"){
      return(bb)
    } else if(type == "gamma"){
      return(gg)
    }

  }

}
