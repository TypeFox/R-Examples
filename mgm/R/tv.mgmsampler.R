

tv.mgmsampler <- function(type, # p time vector
                          lev, # p level vector
                          graphs, # p x p x timestep array
                          threshs, # n list with each p threshold entries (same structure as in mgmsampler)
                          parmatrices = NA,
                          nIter = 250, 
                          varadj = .2) {
  
  # edge weight or parameter matrix?
  if(!is.na(parmatrices)) {
    graphs <- NA
    # basic input info
    n <- dim(parmatrices)[3]
    p <- length(type)
  } else {
    # basic input info
    n <- dim(graphs)[3]
    p <- dim(graphs)[2]
    # input checks
    if(sum(!apply(graphs, 3, matrixcalc::is.symmetric.matrix))>0) stop('The weight matrix must be symmetric for every time step.')
    if(length(threshs)!=n) stop('A threshold for each time step has to be specified.')

  }
  
  # sampling 
  data <- matrix(NA, n, p)
  for(timestep in 1:n) {
    
    if(!is.na(parmatrices)) {
      data[timestep,] <- mgmsampler(n=1, type=type, lev=lev, graph=NA, 
                                    thresh = threshs[[timestep]], parmatrix = parmatrices[[timestep]], 
                                    nIter = nIter, varadj = varadj)
    } else {
      data[timestep,] <- mgmsampler(n=1, type=type, lev=lev, graph=graphs[,,timestep], 
                                    thresh = threshs[[timestep]], parmatrix = NA, 
                                    nIter = nIter, varadj = varadj)
    }
    

  }
  
  # output
  return(data)
  
  
}






