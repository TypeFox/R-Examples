gw.zi.cv <- function(bw, formula, family, dframe, obs, kernel, dmatrix)
{  
  Gl.Model <- eval(substitute(zeroinfl(formula, data = dframe, dist=family)))
  
  RNames<-names(Gl.Model$coefficients$count)
  
  CountNames<-paste("CM",RNames, sep="_")
  ZeroNames<-paste("ZM",RNames, sep="_")
  
  ModelVarNo<-length(Gl.Model$coefficients$count)
  
  ZI_LEst_count<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), CountNames[1:ModelVarNo]))
  ZI_LEst_zero<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), ZeroNames[1:ModelVarNo]))
  ZI_GofFit<-data.frame(y=numeric(0), ZI_yfit=numeric(0), ZI_Res=numeric(0))

  if(kernel == 'adaptive')
    {
       Ne <- bw
    } 
  
  for(m in 1:obs){
      
      #Get the data
      DNeighbour <- dmatrix[,m]
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
      
      #Leave-one-out
      Wts[1]=0
      
      #Calculate WGLM
      ZI.Model <- eval(substitute(zeroinfl(formula, data = SubSet, weights=Wts, dist=family)))

      #Store in table
      ZI_LEst_count[m,]<-ZI.Model$coefficients$count
      ZI_LEst_zero[m,]<-ZI.Model$coefficients$zero

      ZI_GofFit[m,1]<-Gl.Model$y[m]
      ZI_GofFit[m,2]<-exp(sum(Gl.Model$model[m,2:ModelVarNo] * ZI.Model$coefficients$count[2:ModelVarNo]) + sum(ZI.Model$coefficients$count[1]))
   }
  
  Adj.f.l <- sum(ZI_GofFit[,1])/sum(ZI_GofFit[,2])
  ZI_GofFit[,3] <- ZI_GofFit[,1] - Adj.f.l*ZI_GofFit[,2]
  
  RSS <- sum(ZI_GofFit$ZI_Res^2)

  return(RSS)
}

gw.zi.bw<-function(formula, family, dframe, coords, kernel, algorithm="exhaustive", optim.method="Nelder-Mead", b.min=NULL, b.max=NULL, step=NULL) {
  
  Obs <- nrow(dframe)
  DistanceTable <- dist(coords)
  Dij <- as.matrix(DistanceTable)
    
  if(kernel == 'adaptive') 
    {
      if (is.null(b.max) || b.max > Obs) {b.max <- Obs}
      if (is.null(b.min)) {b.min <- 30} 
      
      b <- b.min:b.max
    }
  else 
    { 
      if(kernel == 'fixed')
        {
        Dij.max <- max(Dij)
        Dij.min <- min(Dij[Dij>0])
        
        if (is.null(b.max) || b.max > Dij.max) {b.max <- Dij.max}
        if (is.null(b.min)) {b.min <- Dij.min}  
        if (is.null(step)) {step <- (b.max - b.min)/99}
        
        b<-seq(from=b.min, to=b.max+step, by=step)
        }
    } 
  cat("\nNumber of Observations:", Obs)

  if (algorithm == "exhaustive")
    {
        CVs <- matrix(data=NA, nrow=length(b), ncol=2)
        
        counter <- 1
        
        for(bw in b)
          {
          CV <- try(gw.zi.cv(bw, formula, family, dframe, Obs, kernel, Dij))
          
          if(inherits(CV, "try-error"))
          {
            cat("\nAn error occurred, iteration skipped for bandwidth:", bw)
            CVs[counter,1] <- bw
            CVs[counter,2] <- NA
            counter<-counter+1
            next
          }
          cat("\nBandwidth: ", bw, "CV: ", CV)
          CVs[counter,1] <- bw
          CVs[counter,2] <- CV
          counter<-counter+1
        }
        
        plot(CVs[,1],CVs[,2], xlab="Bandwidth", ylab="Cross-Validation score")
        CV<-min(CVs[,2], na.rm = TRUE)
        bw <- round(CVs[which(CVs[,2]==CV),1],0)
        out<-list(bw=bw, CV=CV, CVs=CVs)
      } 
  else 
  { 
    if (algorithm == "optim") 
      { if(optim.method=="Brent")
      
        {
        opt <- optim(b.min, gw.zi.cv, formula=formula, family=family, dframe=dframe, kernel=kernel, obs=Obs,
                     dmatrix=Dij, method=optim.method, lower = b.min, upper = b.max)        
        
        }
        else{
          opt <- optim(b.min, gw.zi.cv, formula=formula, family=family, dframe=dframe, kernel=kernel, obs=Obs,
                       dmatrix=Dij, method=optim.method)        
      }

        
      }
      bw <- round(opt$par,0)
     out <- list(bw=bw, CV=opt$value)
    #out<-opt
  }

  return(out)
}
