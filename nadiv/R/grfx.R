grfx <- function(n, G, incidence = NULL, saveIncidence = FALSE, output = "matrix", stdnorms = NULL, warn = TRUE){
  d <- nrow(G)
  if(d > 1 && all(G == G[1,1])) warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
  Mg <- as(chol(G), "dtCMatrix")
  if(is.null(incidence)){
     if(any(ls(envir = globalenv() ) == "nadiv_prev_Mincidence")){
       if(warn) warning("using previous incidence matrix")
     } else{
          if(saveIncidence){
	     nadiv_prev_Mincidence <<- Diagonal(n, 1)
             } else{
                  nadiv_prev_Mincidence <- Diagonal(n, 1)
               }
          if(warn) warning("Incidence matrix used = Identity matrix")
       }
  } else{
       if(saveIncidence){
          nadiv_prev_Mincidence <<- chol(incidence)
       } else{
            nadiv_prev_Mincidence <- chol(incidence)
         }
    }

  M <- suppressMessages(kronecker(nadiv_prev_Mincidence, Mg))
  if(is.null(stdnorms)){
     Z <- Matrix(rnorm(n*d), nrow = 1)
  } else{
       if(length(stdnorms) != n*d) stop("length(stdnorms) must be equal to 'n' times the order of 'G'")
       Z <- Matrix(stdnorms, nrow = 1)
    }
  X <- Matrix((Z %*% M)@x, ncol = d, byrow = TRUE)

 return(as(X, output))
}

