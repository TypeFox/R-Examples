FUNTA <-
function(Data, centered = FALSE, give.angles = FALSE, tick.dist = 1){
  # Data: a matrix with n rows (observations) and p columns (time points)
  # centered: Boolean. If FALSE, data are centered before calculating intersection angles.
  #                    If TRUE, data are assumed to be centered. Warning: If data are not centered
  #                    and TRUE is chosen, the computation fails if one curve has no intersections.
  # give.angles: Boolean. If FALSE, the intersection angles are not shown. If TRUE, they are shown.
  # tick.dist: numeric. The distance between two time points data can be specified here. 
  nc <- ncol(Data)
  nr <- nrow(Data)
  ftad <- numeric(nr)
  
  if (centered == FALSE) Data <- Data - rowMeans(Data) # quick and dirty centering of Data
  if (give.angles == TRUE) angle.near2 <- list()
  
  for (i in 1:nr){ # for each of the curves, do
    DataC <- matrix(Data[-i,], ncol = ncol(Data))
    # 1. Find intersections
    # intersections on given time points
    if (any(Data[i,] == t(DataC)) == TRUE){ # if at least one intersection on 
      # given time points exist, 
      loc.exact <- which(Data[i,] == t(DataC), arr.ind=TRUE) # save the location
      count.exact <- nrow(loc.exact)
      angle.exact <- numeric(count.exact)
      
      # 2. Calculate the angles
      for (j in 1:count.exact){
        if (loc.exact[j,1] < nc & loc.exact[j,1] > 1){
          cut.slope <- diff(DataC[loc.exact[j,2], 
                                  c((loc.exact[j,1] - 1), (loc.exact[j,1] + 1)) ])
          ref.slope <- diff(Data[i, c((loc.exact[j,1] - 1), (loc.exact[j,1] + 1)) ])
          angle.exact[j] <- angle(cut.slope, ref.slope, tick.dist)
        }
        else angle.exact[j] <- 0
      }
    }
    
    # 1. Find intersections    
    # Intersections near given time points
    sign.near <- sign(t(Data[i,] - t(DataC)))
    loc.near <- which(sign.near[,-nc] - sign.near[,-1] != 0, arr.ind = TRUE) 
    # switches from minus to plus = -2, 
    # from plus to minus = 2
    if (!is.matrix(loc.near)) loc.near <- cbind(1, loc.near) # 1 because ow. the index would be missing
    
    bb <- nrow(loc.near)
    angle.near <- numeric(bb) 
    
    # 2. Calculate the angles   
    for (j in 1:bb){    
      cut.slope <- diff(DataC[loc.near[j,1], (loc.near[j,2] : (loc.near[j,2] + 1))])
      ref.slope <- diff(Data[i, (loc.near[j,2] : (loc.near[j,2] + 1))])
      #      if (give.angles == FALSE) 
      angle.near[j] <- angle(cut.slope, ref.slope, tick.dist)
      #      else angle.near2[j] <- angle(cut.slope, ref.slope, tick.dist)
    }
    if (give.angles == TRUE) angle.near2[[i]] <- angle.near
    
    # combine results into a single vector
    if (any(Data[i,] == t(DataC)) == TRUE){ # if at least one intersection on 
      if (give.angles == FALSE) angle.near <- c(angle.exact, angle.near)
      else angle.near2[[i]] <- c(angle.exact, angle.near)
      #        angle.exact <- numeric()
    }
    
    # 3. Calculate the depth
    if (give.angles == FALSE) ftad[i] <- 1 - mean(angle.near)/pi
    else ftad[i] <- 1 - lapply(FUN = mean, angle.near2)[[i]]/pi
  }
  if (give.angles == FALSE) ftad
  else list(FUNTA = ftad, Angles = lapply(angle.near2, "*", 180/pi))
}
