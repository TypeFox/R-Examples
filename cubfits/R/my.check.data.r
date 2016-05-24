### This function will check if input data are all well-behaved.

my.check.data <- function(phi.Obs = NULL, phi.Init = NULL,
    phi.pred.Init = NULL){
  ### Check phi.* are well-behaved.
  if(!is.null(phi.Obs)){
    if(!(all(is.finite(phi.Obs)) && all(phi.Obs > 0))){
      .cubfitsEnv$my.stop("phi.Obs is invalid.")
    }
    if(.CF.CONF$scale.phi.Obs){
      if(abs(mean(phi.Obs) - 1) > 1e-8){
        .cubfitsEnv$my.stop(paste("mean(phi.Obs) =", mean(phi.Obs)))
      }
    }
  }
  if(!is.null(phi.Init)){
    if(!(all(is.finite(phi.Init)) && all(phi.Init > 0))){
      .cubfitsEnv$my.stop("phi.Init is invalid.")
    }
    if(.CF.CONF$scale.phi.Obs || .CF.CONF$estimate.bias.Phi){
      if(abs(mean(phi.Init) - 1) > 1e-8){
        .cubfitsEnv$my.stop(paste("mean(phi.Init) =", mean(phi.Init)))
      }
    }
  }
  if(!is.null(phi.pred.Init)){
    if(!(all(is.finite(phi.pred.Init)) && all(phi.pred.Init > 0))){
      .cubfitsEnv$my.stop("phi.pred.Init is invalid.")
    }
    if(.CF.CONF$scale.phi.Obs || .CF.CONF$estimate.bias.Phi){
      if(abs(mean(phi.pred.Init) - 1) > 1e-8){
        .cubfitsEnv$my.stop(paste("mean(phi.pred.Init) =", mean(phi.pred.Init)))
      }
    }
  }

  invisible()
} # End of my.check.data().
