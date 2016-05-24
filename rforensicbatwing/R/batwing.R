#reps
#Number of output lines

#burnin
#Number of reps to take before starting recording data

#Nbetsamp
#The number of times that changes to hyperparameters are attempted between outputs.

#treebetN 
#The number of times that changes to the genealogical tree are attempted
#before any changes to the hyperparameters are attempted. Thus BATWING
#outputs are separated by treebetN x Nbetsamp attempted tree updates.

#sizemodel
#Code for the population growth model: 
#0, constant population size; 
#1, exponential growth at rate alpha at all times;
# alphaprior = NULL gives sizemodel = 0, i.e. constant population size
# alphaprior a valid string gives exponential population size

batwing <-
function(database, reps = 10, burnin = 0, treebetN = 10, Nbetsamp = 10, muprior = "constant(0.003)", Nprior = "lognormal(9, 1)", alphaprior = NULL, progress = TRUE, trace = FALSE) { 

  if (any(database <= 0)) {
    stop("database must consist of positive integers (> 0) only")
  }
  
  res <- .internalbatwing(database = database, 
    reps = reps, 
    burnin = burnin, 
    treebetN = treebetN, 
    Nbetsamp = Nbetsamp, 
    muprior = muprior, 
    Nprior = Nprior, 
    alphaprior = alphaprior, 
    forensicmode = FALSE, 
    progress = progress, trace = trace)
  
  return(res)
}

