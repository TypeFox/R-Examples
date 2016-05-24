eigenvalueBased <- function(g, matrix_function, s=1){  
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))

  M <- do.call(matrix_function,list(g))

  EV <- as.double(abs(eigen(M,only.values=TRUE)$values));
  EV <- EV[EV != 0];
  EVs <- EV^(1/s)
  sumEVs <- sum(EVs)
  pi<- EVs/sumEVs

  result <- list()

  ##Expression (2)
  result[["HMs"]] <- (-1) * sum(pi*log2(pi))

  ##Expression (3)
  result[["SMs"]] <- sumEVs

  ##Expression (4)
  result[["ISMs"]] <- 1/sumEVs

  ##Expression (5)
  result[["PMs"]] <- prod(EVs)

  ##Expression (6)
  result[["IPMs"]] <- 1/result[["PMs"]]

  result
}
