Update_Lambda_Estimates <- function(i,
                                    gpar,
                                    theta,
                                    together,
                                    verbose,
                                    net,
                                    GERGM_Object){
  # if we are usinga correlation network, do beta regression
  if(GERGM_Object@is_correlation_network){
    stop("Currently not implemented! Set omit_intercept_term = TRUE and include an edges term in specification...")
  }else{
    #do our normal t regression
    cat("Updating Estimates -- Iteration:", i," \n")
    GERGM_Object <- store_console_output(GERGM_Object,paste("Updating Estimates -- Iteration:", i," \n"))
    if(verbose){
      cat("Lambda Estimates", gpar$par,"\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object,paste("Lambda Estimates", gpar$par,"\n"))
    gpar.new <- NULL
    if(verbose){
      gpar.new <- optim(par = as.numeric(gpar$par),
                        llg,
                        alpha = GERGM_Object@weights,
                        theta = as.numeric(theta$par),
                        z = GERGM_Object@data_transformation,
                        method = "BFGS",
                        together = together,
                        GERGM_Object = GERGM_Object,
                        hessian = T,
                        control = list(fnscale = -1, trace = 6))
    }else{
      gpar.new <- optim(par = as.numeric(gpar$par),
                        llg,
                        alpha = GERGM_Object@weights,
                        theta = as.numeric(theta$par),
                        z = GERGM_Object@data_transformation,
                        method = "BFGS",
                        together = together,
                        GERGM_Object = GERGM_Object,
                        hessian = T,
                        control = list(fnscale = -1, trace = 0))
    }
    if(verbose){
      cat("Lambda estimates", "\n")
    }
    GERGM_Object <- store_console_output(GERGM_Object, "Lambda estimates\n")
    if(verbose){
      print(gpar.new$par)
    }
    GERGM_Object <- store_console_output(GERGM_Object, toString(gpar.new$par))
    gpar.std.errors <- 1 / sqrt(abs(diag(gpar.new$hessian)))
    # Transform the unbounded weights to bounded weights via a t-distribution
    beta <- gpar.new$par[1:(length(gpar.new$par) - 1)]
    sig <- 0.01 + exp(gpar.new$par[length(gpar.new$par)])
    BZ <- 0
    for (j in 1:(dim(GERGM_Object@data_transformation)[3])) {
      BZ <- BZ + beta[j] * GERGM_Object@data_transformation[, , j]
    }

    #store so we can transform back
    GERGM_Object@BZ <- BZ
    GERGM_Object@BZstdev <- sig
    transformation_type <- GERGM_Object@transformation_type
    if(transformation_type == "logcauchy"){
      GERGM_Object@bounded.network <- pst(log(net), BZ, sig, 1)
    }
    if( transformation_type == "cauchy"){
      GERGM_Object@bounded.network <- pst(net, BZ, sig, 1)
    }
    if(transformation_type == "lognormal"){
      GERGM_Object@bounded.network <- pst(log(net), BZ, sig, Inf)
    }
    if( transformation_type == "gaussian"){
      GERGM_Object@bounded.network <- pst(net, BZ, sig, Inf)
    }
  }
  return(list(GERGM_Object = GERGM_Object,
              gpar.new = gpar.new,
              gpar.std.errors  = gpar.std.errors))
}


