gw.glm.cv <- function(bw, formula, family, dframe, obs, kernel, dmatrix)
{  
  Gl.Model <- eval(substitute(glm(formula, family=family, data = dframe)))
  RNames<-names(Gl.Model$coefficients)
  ModelVarNo <-length(Gl.Model$coefficients)
  
  GGLM_LEst<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), RNames[1:ModelVarNo]))
  GGLM_GofFit<-data.frame(y=numeric(0), GLM_yfit=numeric(0), GLM_Res=numeric(0))
  
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
      Lcl.Model<-eval(substitute(glm(formula, family=family, data = SubSet, weights=Wts)))

      #Store in table
      GGLM_LEst[m,]<-Lcl.Model$coefficients
      GGLM_GofFit[m,1]<-Gl.Model$y[m]
      GGLM_GofFit[m,2]<-exp(sum(Gl.Model$model[m,2:ModelVarNo] * Lcl.Model$coefficients[2:ModelVarNo]) + sum(Lcl.Model$coefficients[1]))
    }
  
  Adj.f.l <- sum(GGLM_GofFit[,1])/sum(GGLM_GofFit[,2])
  GGLM_GofFit[,3] <- GGLM_GofFit[,1] - Adj.f.l*GGLM_GofFit[,2]
  
  RSS <- sum(GGLM_GofFit$GLM_Res^2)

  return(RSS)
}

gw.glm.bw<-function(formula, family, dframe, coords, kernel, algorithm="exhaustive", optim.method="Nelder-Mead", b.min=NULL, b.max=NULL, step=NULL) {
  
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
          CV <- gw.glm.cv(bw, formula, family, dframe, Obs, kernel, Dij)
          cat("\nBandwidth: ", bw, "CV: ", CV)
          CVs[counter,1] <- bw
          CVs[counter,2] <- CV
          counter<-counter+1
          }
        plot(CVs[,1],CVs[,2], xlab="Bandwidth", ylab="Cross-Validation score")
        CV<-min(CVs[,2])
        bw <- round(CVs[which(CVs[,2]==CV),1],0)
        out<-list(bw=bw, CV=CV, CVs=CVs)
      } 
  else 
  { 
    if (algorithm == "optim") 
      { if(optim.method=="Brent")
      
        {
        opt <- optim(b.min, gw.glm.cv, formula=formula, family=family, dframe=dframe, kernel=kernel, obs=Obs,
                     dmatrix=Dij, method=optim.method, lower = b.min, upper = b.max)        
        
        }
        else{
          opt <- optim(b.min, gw.glm.cv, formula=formula, family=family, dframe=dframe, kernel=kernel, obs=Obs,
                       dmatrix=Dij, method=optim.method)        
      }

        
      }
      bw <- round(opt$par,0)
     out <- list(bw=bw, CV=opt$value)
    #out<-opt
  }

  return(out)
}
