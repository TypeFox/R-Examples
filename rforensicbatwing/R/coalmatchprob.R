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

coalmatchprob <-
function(database, haplotype, reps = 10, burnin = 0, treebetN = 10, Nbetsamp = 10, muprior = "constant(0.003)", Nprior = "lognormal(9, 1)", alphaprior = NULL, progress = TRUE, trace = FALSE) { 
  if (is.null(database) || !is.matrix(database) || ncol(database) == 0 || nrow(database) == 0 || !is.numeric(database)) {
    stop("database must be a matrix of integers")
  }
  if (nrow(database) == 1) {
    database <- matrix(as.integer(database), nrow = 1)
  } else {
    database <- apply(database, 2, as.integer)
  }
  
  if (is.null(haplotype) || !is.numeric(haplotype) || length(haplotype) != ncol(database)) {
    stop("haplotype must have a length corresponding to the number of columns of database")
  }
  haplotype <- as.integer(haplotype)  


  if (any(database <= 0)) {
    stop("database must consist of positive integers (> 0) only")
  }
  if (any(haplotype <= 0)) {
    stop("haplotype must consist of positive integers (> 0) only")
  }


  missing <- rep(-1, length(haplotype))
  data <- rbind(haplotype, missing, database)
  attr(data, "dimnames") <- NULL
  
  res <- .internalbatwing(database = data, 
    reps = reps, 
    burnin = burnin, 
    treebetN = treebetN, 
    Nbetsamp = Nbetsamp, 
    muprior = muprior, 
    Nprior = Nprior, 
    alphaprior = alphaprior, 
    forensicmode = TRUE, 
    progress = progress, trace = trace)
  
  return(res)
}

