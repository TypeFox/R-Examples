makedesignkey <- function(keys, primes){
## A function to create an object of class designkey directly
## from a list of simple matrices
  ## ARGUMENTS
  ## - keys : a list of n (integer) matrices
  ## - primes : a vector of n prime numbers
  ## OUTPUT
  ## an object of class designkey
  ## DETAILS
  ## the names of the factors are extracted from the matrix column names
  ## the associated model formula is the additive model with all factors
  ## the estimate formula is the constant
  ## EXAMPLE
  ## mat1 <- cbind(diag(3),1)
  ## colnames(mat1) <- c("A","B","C","D")
  ## mat2 <- cbind(diag(2),c(1,2))
  ## colnames(mat2) <- c("E","F","G")
  ## mat.dk <- makedesignkey(list(mat1), primes=c(2))
  ## print(mat.dk)
  ## summary(mat.dk)
  ## alias(mat.dk)
  ## mat.plan <- planor.design(mat.dk)
  ## -----------------------------------------------------
  ## information extraction
  fact.names <- unlist( lapply(keys, colnames) )
  fact.levels <- rep( primes, unlist(lapply(keys, ncol)) )
  tp.fact <- planor.factors(factors=fact.names, nlevels=fact.levels)
  tp.fact@fact.info$model <- TRUE
  tp.fact@pseudo.info$model <- TRUE
  ##
  nunits <- prod( primes^unlist(lapply(keys, nrow)) )
  ##
  names(keys) <- as.character(primes)
  for(p in seq(keys)){
    keys[[p]] <- new("keymatrix", keys[[p]], p=primes[p])
  }
  ##
  tp.model <- as.formula(paste("~", paste(fact.names,collapse="+"), sep=""))
  tp.estim <- as.formula(paste("~", paste(fact.names,collapse="+"), sep=""))
  tp.mod <- planor.model(model=tp.model, estimate=tp.estim)

  ## KEY GENERATION
  newkey <- new("designkey",
                .Data=keys,
                factors=tp.fact,
                nunits=nunits,
                model=tp.mod)

  return(newkey)
}
