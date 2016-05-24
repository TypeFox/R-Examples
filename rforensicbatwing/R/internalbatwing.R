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

.internalbatwing <-
function(database, reps = 10, burnin = 0, treebetN = 10, Nbetsamp = 10, muprior = "constant(0.003)", Nprior = "lognormal(9, 1)", alphaprior = NULL, forensicmode, progress = TRUE, trace = FALSE) { 
  if (is.null(database) || !is.matrix(database) || ncol(database) == 0 || nrow(database) == 0 || !is.numeric(database)) {
    stop("database must be a matrix of positive integers (> 0)")
  }
  if (nrow(database) == 1) {
    database <- matrix(as.integer(database), nrow = 1)
  } else {
    database <- apply(database, 2, as.integer)
  }
  
  if (is.null(reps) || length(reps) != 1 || !is.numeric(reps) || reps <= 0) {
    stop("reps must be >= 1")
  }
  reps <- as.integer(reps)
  
  if (is.null(burnin) || length(burnin) != 1 || !is.numeric(burnin) || burnin < 0) {
    stop("burnin must be >= 0")
  }
  burnin <- as.integer(burnin)
  
  if (is.null(treebetN) || length(treebetN) != 1 || !is.numeric(treebetN) || treebetN <= 0) {
    stop("treebetN must be >= 1")
  }
  treebetN <- as.integer(treebetN)

  if (is.null(Nbetsamp) || length(Nbetsamp) != 1 || !is.numeric(Nbetsamp) || Nbetsamp <= 0) {
    stop("Nbetsamp must be >= 1")
  }
  Nbetsamp <- as.integer(Nbetsamp)
  
  if (is.null(muprior) || !is.character(muprior)) {
    stop("Specify at least one mu prior as a character vector e.g. 'muprior = \"normal(0.003, 0.001)\"'")
  }
  
  if (length(muprior) != 1 && length(muprior) != ncol(database)) {
    stop("Specify either one mu prior or as many as there are loci")
  }

  if (is.null(Nprior) || !is.character(Nprior) || length(Nprior) != 1) {
    stop("Specify one N prior as a character vector e.g. 'Nprior = \"normal(20000, 3000)\"'")
  }

  if (!is.null(alphaprior) && (!is.character(alphaprior) || length(alphaprior) != 1)) {
    stop("Specify alpha prior as NULL or as a character vector e.g. 'alphaprior = \"uniform(0, 0.01)\"'")
  }

  if (missing(forensicmode) || is.null(forensicmode) || !is.logical(forensicmode)) {
    stop("forensicmode must be logical/boolean")
  }
  
  if (is.null(progress) || !is.logical(progress)) {
    stop("progress must be logical/boolean")
  }
  
  if (is.null(trace) || !is.logical(trace)) {
    stop("trace must be logical/boolean")
  }
  
  ################################################
  
  attr(database, "dimnames") <- NULL
  
  res <- .Call("batwing", database, reps, burnin, treebetN, Nbetsamp, muprior, length(muprior), Nprior, alphaprior, forensicmode, progress, trace, PACKAGE = "rforensicbatwing")
  
  col.i <- colnames(res$result) == "mu"
  if (any(col.i) && sum(col.i) > 1) {
    colnames(res$result)[col.i] <- paste(colnames(res$result)[col.i], 1:sum(col.i), sep = "")
  }
  
  res$result <- data.frame(res$result)
  
  if (forensicmode) {
    class(res) <- c("forensicbatwing", "batwing")
  } else {
    res$result$p <- NULL
    class(res) <- "batwing"
  }
  
  return(res)
}

