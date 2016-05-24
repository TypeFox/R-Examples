############################################################################
# Generate data.frame with simulated data for multinomial logit model
# Args: 
#   modelType <- 'X'    # Set to ONE of 'X', 'Y', 'Z' 
#   K         <- 100    # num choices
#   numCovars <- 50     # number of covariates   
#   obsFactor <- 15     # num obs = K * numCovars * obsFactor
#   outfile   <- FALSE  # Whether to write data.frame object to a .rds file
#   print.level <- 0    # Set to non-zero to print stuff
#   seed        <- pi   # Value of random seed
# Returns:
#   data.frame in 'long' format such that data for a choice is together
############################################################################
makeModel <-function(modelType, K=100, numCovars=50, obsFactor=20,
                     outfile=FALSE, print.level=0, seed=pi) 
{
  set.seed(seed)

  # Automatic params
  stopifnot(sum(modelType == c('X', 'Y', 'Z')) > 0)
  N         <- numCovars*K*obsFactor   # num observations
  outfile <- if (!outfile) NULL 
             else paste0("mlogitSim_", modelType, "_K", K, ".rds") # Filename

  if (print.level) {
      cat(paste("\nGenerating multinomial logit model with:", numCovars, 
                "variables \tN =", N, "observations\n\tK =", K, " choices.\n"),
          file=stdout())
  }

  t0 <- proc.time()[3]/60
  # Make data matrix & coefficients vector & compute probabilities
  sim.df <- if (modelType == 'X') 
      matrix(runif(N*numCovars), nrow=N*K, ncol=numCovars, byrow=TRUE)
  else 
      matrix(runif(N*K*numCovars), nrow=N*K, ncol=numCovars)

  # Generate response
  each <- floor(N/(K*2))
  stopifnot(each >= 1)
  rest <- N - each * K
  ff <- function(i) {
      x <- rep(0, K)
      x[i] <- 1
      rep(x, each) # Generate same response for 'each' number of observations
  }
  gg <- function(i) {
      x <- rep(0, K)
      x[sample(1:K, 1)] <- 1
      x 
  }
  respVec <- c(as.vector(sapply(seq(1, K), ff)),
               as.vector(sapply(seq(1, rest), gg)))
  mincount <- min(colSums(matrix(respVec, nrow=N, ncol=K)))

  # Generate data.frame. Insert: response, choices, indivID  
  sim.df <- data.frame(sim.df)
  formulaStr <- paste0("response ~ 0 + ", paste(colnames(sim.df), collapse=" + ")) 
  choiceNames <- rep(paste0("ch", seq(1, K)), N)
  ordering <- order(choiceNames)
  sim.df$choices <- sort(choiceNames)
  sim.df$response <- respVec[ordering]  
  sim.df$indivID <- rep(seq(1, N), K)
     
  t1 <- proc.time()[3]/60
  if (print.level)
      cat(paste("Time taken for making data.frame (min) ="), t1 - t0, "\n")
  
  # Save to file
  if (!is.null(outfile)) {
    cat(paste0("Writing data.frame object: 'sim.df' to file: "), outfile, "\n")
    saveRDS(sim.df, file=outfile)
  }
  if (print.level) {
      cat(paste0("Lowest frequency of any choice is = ", mincount/N, "\n"))
      cat(paste0("Model formula is: \n", formulaStr, "\n")) 
  }
  sim.df
}
