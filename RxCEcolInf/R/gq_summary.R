`gq.summary` <- function(gq.output,burnin=0,quantiles=c(0.01,0.025,0.25,0.50,0.75,0.975,0.99))
{
  # Convert the GQ output to an mcmc.list:
  if (!is.null(gq.output$draws)){
    # Should be a single gq.output:
    ndraws <- nrow(gq.output$draws$mu)
    gq.output <- list(gq.output)
  } else {
    # gq.draws should be a list of gq.outputs:
    if (is.null(gq.output[[1]]$draws)){
      stop("gq.output must be either a list of or a single return object from 'Analyze'")
    }
    ndraws <- nrow(gq.output[[1]]$draws$mu)
  }
  if ((burnin<0)||(burnin>=ndraws)){
    stop("Invalid burn-in period")
  }
  # Create the mcmc.list (overwrite first element):
  gq.draws <- mcmc.list(mcmc(1))

  for (i in 1:length(gq.output)){

    # Create a matrix of draws from the i^th chain:
    tmp.chain <- cbind(gq.output[[i]]$draws$mu[(burnin+1):ndraws,],
		       gq.output[[i]]$draws$Sigma[(burnin+1):ndraws,],
		       gq.output[[i]]$draws$NNs[(burnin+1):ndraws,],
		       gq.output[[i]]$draws$LAMBDA[(burnin+1):ndraws,],
		       gq.output[[i]]$draws$TURNOUT[(burnin+1):ndraws,],
		       gq.output[[i]]$draws$GAMMA[(burnin+1):ndraws,])

    # Add that chain to the mcmc.list:
    gq.draws[[i]] <- mcmc(tmp.chain)
  }
  return(summary(gq.draws,quantiles=quantiles))
}
