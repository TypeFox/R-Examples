gw.zi.light <- function(formula, family, dframe, bw, kernel, coords)
{  
  Obs <- nrow(dframe)
  if(kernel == 'adaptive')
  { 
    Ne <- bw
    #cat("\nKernel: Adaptive\nNeightbours:", Ne)
  } 
  else 
  { 
    if(kernel == 'fixed')
    {
     # cat("\nKernel: Fixed\nBandwidth:", bw)
    }
  }
  
  Gl.Model <- eval(substitute(zeroinfl(formula, data = dframe, dist=family)))
  
  RNames <- names(Gl.Model$coefficients$count)
  CountNames <- paste("CM",RNames, sep="_")
  ZeroNames <- paste("ZM",RNames, sep="_")

  ModelVarNo <- length(Gl.Model$coefficients$count)
  
  DistanceT <- dist(coords)
  Dij <- as.matrix(DistanceT)
  
  ZI_LEst_count<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), CountNames[1:ModelVarNo]))
  ZI_LEst_zero<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), ZeroNames[1:ModelVarNo]))

  for(m in 1:Obs){
    
    #Get the data
    DNeighbour <- Dij[,m]
    DataSet <- data.frame(dframe, DNeighbour=DNeighbour)
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    if(kernel == 'adaptive')
    { 
      #Keep Nearest Neighbours
      SubSet <- DataSetSorted[1:Ne,]
      Kernel_H <- max(SubSet$DNeighbour)
    } 
    else 
    { 
      if(kernel == 'fixed')
      {
        SubSet <- subset(DataSetSorted, DNeighbour <= bw)
        Kernel_H <- bw
      }
    }
    
    #Bi-square weights
    Wts<-(1-(SubSet$DNeighbour/Kernel_H)^2)^2
    
    #Calculate GW-ZI
    ZI.Model <- eval(substitute(zeroinfl(formula, data = SubSet, weights=Wts, dist=family)))
    
    #print(summary(ZI.Model))
    
    #Store in table
    
    ZI_LEst_count[m,]<-ZI.Model$coefficients$count
    ZI_LEst_zero[m,]<-ZI.Model$coefficients$zero
  }
  
  results<-list(ZI_LEst_count=ZI_LEst_count, ZI_LEst_zero=ZI_LEst_zero)
  return(results)
}

#-----------Monte Carlo Test --------------------------
gw.zi.mc.test <- function (Nsim = 19, formula, family, dframe, bw, kernel, coords) 
{
  
  Obs <- nrow(dframe)
  
  if (family=="poisson") 
  {
    
    gw.zi.observed <- eval(substitute(gw.zi.light(formula, family, dframe, bw, kernel, coords)))
    
    var.lpest.count.obs <- diag(var(gw.zi.observed$ZI_LEst_count))
    var.lpest.zero.obs <- diag(var(gw.zi.observed$ZI_LEst_zero))
    
#     print("Observed")
#     print(round(var.lpest.count.obs,5)) 
#     print(round(var.lpest.zero.obs,5))
    
    params<-length(var.lpest.count.obs)
    
    l.pest.SIM <- matrix(data = NA, nrow = Nsim, ncol = 2*params+1)
    l.pest.SIM.c<- matrix(data = 0, nrow = Nsim, ncol = 2*params+1)
  }
  
  for (i in 1:Nsim) {
    x = runif(Obs, min = 1, max = Obs * (Nsim + 1))
    tx <- round(x/(Nsim + 1))
    CoordsTmp <- data.frame(coords, RandomID = tx)
    CoordsSorted <- CoordsTmp[order(CoordsTmp$RandomID),]
    CoordsNew <- CoordsSorted[,1:2]
    
    gw.zi.SIM <- try(eval(substitute(gw.zi.light(formula, family, dframe, bw, kernel, CoordsNew))))
    
    if(inherits( gw.zi.SIM, "try-error"))
    {
      #on error repeat this iteration using
      #NAs
      i <- i-1
      cat("\nAn error occurred, iteration to be repeated:", i)
      next
      
      
    }
    
    l.pest.SIM[i, 1] <- i
    l.pest.SIM.c[i, 1] <- i
    
    for (j in 1:params) {
      
      l.pest.SIM[i, j+1] <- var(gw.zi.SIM$ZI_LEst_count[,j])
      
      if (l.pest.SIM[i, j+1] >= var.lpest.count.obs[[j]]) {l.pest.SIM.c[i, j+1] <- 1}
    }
    
    for (j in (params+1):(2*params)) {
      
      l.pest.SIM[i, j+1] <- var(gw.zi.SIM$ZI_LEst_zero[,j-params])
      
      if (l.pest.SIM[i, j+1] >= var.lpest.zero.obs[[j-params]]) {l.pest.SIM.c[i, j+1] <- 1}
    }
  }
  
  C.test <- colSums(l.pest.SIM.c)
  
  for (j in 1:(2*params)) {
    if ((Nsim - C.test[[j+1]]) < C.test[[j+1]]) {
      C.test[[j+1]] = Nsim - C.test[[j+1]]
    }
  }
  
  pseudo.p = (1 + C.test)/(Nsim + 1)
  
  return(list(var.lpest.obs = cbind(var.lpest.count.obs,var.lpest.zero.obs), var.SIM = l.pest.SIM, var.SIM.c = l.pest.SIM.c, pseudo.p = pseudo.p[2:(2*params+1)]))
  
}

