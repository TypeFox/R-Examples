rFUNTA <-
function(Data, centered = FALSE, type.inner = "max", 
                   type.outer = "median", tick.dist = 1, nObs = nrow(Data)){
  p <- nrow(Data) / nObs
  if(p != round(p)) stop("Error: Wrong number of observations. Choose nObs such that nrow(Data)/nObs is an integer.")
  FUNTA.curr <- matrix(nrow = p, ncol = nObs)
  temp.mi <- matrix(nrow = p, ncol = nObs)
  # partitioning the data
  lower.bounds <- seq(1, nrow(Data), by = nObs)
  upper.bounds <- seq(nObs, nrow(Data), length = p)
  nc <- ncol(Data)
  #  nr <- nrow(Data)
  ftad <- numeric(nObs)
  angle.each <- matrix(nrow = p, ncol = nObs-1) # to store max/mean/.. angle of curve i with curve k in
  
  for (i in 1:nObs){ # for each of the curves, do
    for (j in 1:p){ # for each dimension, do
      curr.dimension <- Data[lower.bounds[j]:upper.bounds[j], ]
      if (centered == FALSE) curr.dimension <- curr.dimension - rowMeans(curr.dimension) # quick and dirty centering of Data 
      
      ind.without.i <- (1:nObs)[-i]
      m <- 1
      
      for (k in ind.without.i){ # look every other curve
        DataC <- curr.dimension[k,]
        # 1. Find intersections
        # intersections on given time points
        if (any(curr.dimension[i,] == DataC) == TRUE){ # if at least one intersection on 
          # given time points exist, 
          loc.exact <- which(curr.dimension[i,] == DataC) # save the location
          count.exact <- length(loc.exact)
          
          angle.exact <- numeric(count.exact)
          
          # 2. Calculate the angles
          for (l in 1:count.exact)
          {
            if (loc.exact[l] < nc & loc.exact[l] > 1){
              cut.slope <- diff(DataC[c((loc.exact[l] - 1), (loc.exact[l] + 1)) ])
              ref.slope <- diff(curr.dimension[l, c((loc.exact[l] - 1), (loc.exact[l] + 1)) ])
              
              angle.exact[l] <- angle(cut.slope, ref.slope, tick.dist)
            }
            else angle.exact[l] <- 0
          }
        }
        
        # 1. Find intersections    
        # Intersections near given time points
        sign.near <- sign(curr.dimension[i,] - DataC)
        
        loc.near <- which(sign.near[-nc] - sign.near[-1] != 0) 
        # switches from minus to plus = -2, 
        # from plus to minus = 2
        
        angle.near <- numeric(length(loc.near))
        
        # 2. Calculate the angles   
        for (l in 1:length(loc.near)){    
          cut.slope <- diff(DataC[loc.near[l] : (loc.near[l] + 1)])
          ref.slope <- diff(curr.dimension[i, (loc.near[l] : (loc.near[l] + 1))])
          
          angle.near[l] <- angle(cut.slope, ref.slope, tick.dist)
        }
        
        # combine results into a single vector
        if (any(curr.dimension[i,] == DataC) == TRUE) angle.near <- c(angle.exact, angle.near) 
        # if at least one intersection on 
        
        
        if (type.inner == "median") angle.each[j,m] <- median(angle.near)
        if (type.inner == "mean") angle.each[j,m] <- mean(angle.near)
        if (type.inner == "max") angle.each[j,m] <- max(angle.near)
        m <- m+1
      }
    }      
    angle.dim.max <- apply(angle.each, 2, max)
    # 3. Calculate the depth
    if (type.outer == "median") ftad[i] <- 1 - median(angle.dim.max) / pi
    if (type.outer == "mean") ftad[i] <- 1 - mean(angle.dim.max) / pi
    if (type.outer == "max") ftad[i] <- 1 - max(angle.dim.max) / pi
    angle.each <- matrix(nrow = p, ncol = nObs-1)
  }
  ftad
}
