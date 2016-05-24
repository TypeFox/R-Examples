
mcmc2bugs <- function(x, model.file = NULL, program = "", DIC = FALSE, 
  DICOutput = NULL, n.iter = NULL, n.burnin = 0, n.thin = 1){
  parameter.names <- dimnames(x[[1]])[[2]]
  n.keeps <- dim(x[[1]])[1]
  n.chains <- summary(x)[["nchain"]]
  n.parameters.to.save <- length(parameter.names)
  sims.array <- array(NA, c(n.keeps, n.chains, n.parameters.to.save))
  dimnames(sims.array)[[3]] <- parameter.names
  for (i in 1:n.chains){
    sims.array[,i,] <- x[[i]]
  }
  ans <- as.bugs.array2(sims.array, model.file=model.file, program=program, DIC=DIC, 
    DICOutput=DICOutput, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin)
  return(ans)
}
