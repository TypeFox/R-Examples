gwr <- function(formula, dframe, bw, kernel, coords)
{  
  Obs <- nrow(dframe)
  
  cat("\nNumber of Observations:", Obs)
  
  if(kernel == 'adaptive')
  { 
    Ne <- bw
    cat("\nKernel: Adaptive\nNeightbours:", Ne)
  } 
  else 
  { 
    if(kernel == 'fixed')
    {
      cat("\nKernel: Fixed\nBandwidth:", bw)
    }
  }
#  VarNo<-ncol(DFrame)
  
  Gl.Model<-eval(substitute(lm(formula, data = dframe)))
  
  RNames<-names(Gl.Model$coefficients)
  PValueNames<-paste("P",RNames, sep="_")
  ModelVarNo<-length(Gl.Model$coefficients)
  
  cat("\nNumber of Variables:", ModelVarNo-1)
  cat("\n--------------- Global Model Summary ---------------\n")
  print(summary(Gl.Model))  
  
  DistanceT <- dist(coords)
  Dij <- as.matrix(DistanceT)
  
  LM_LEst<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), RNames[1:ModelVarNo]))
  LM_LPvalues<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), PValueNames[1:ModelVarNo]))
  LM_GofFit<-data.frame(y=numeric(0), LM_yfit=numeric(0), LM_Res=numeric(0), LM_AIC=numeric(0),LM_Deviance=numeric(0))

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
    
    #Calculate WLM
    Lcl.Model<-eval(substitute(lm(formula, data = SubSet, weights=Wts)))
    
    #Store in table
    LM_LEst[m,]<-Lcl.Model$coefficients
    LM_LPvalues[m,]<-summary(Lcl.Model)$coefficients[,4]
    LM_GofFit[m,1]<-Gl.Model$model[m,1]
    LM_GofFit[m,2]<-sum(Gl.Model$model[m,2:ModelVarNo] * Lcl.Model$coefficients[2:ModelVarNo]) + Lcl.Model$coefficients[[1]]
    LM_GofFit[m,3]<-LM_GofFit[m,1] - LM_GofFit[m,2]
    LM_GofFit[m,4]<-AIC(Lcl.Model)
    LM_GofFit[m,5]<-deviance(Lcl.Model)
  }
  
  gw.lm.out<-list(LM_LEst= LM_LEst, LM_LPvalues=LM_LPvalues, LM_GofFit=LM_GofFit)
  
  cat("\n--------------- Local Model Summary ---------------\n")
 
  cat("\nResiduals:\n")
  print(summary(gw.lm.out$LM_GofFit$LM_Res))
  
  
  t1 <- data.frame(Min = apply(gw.lm.out$LM_LEst, 2, min), Max = apply(gw.lm.out$LM_LEst, 2, max), 
                   Mean = apply(gw.lm.out$LM_LEst, 2, mean), StD = apply(gw.lm.out$LM_LEst, 2, sd))
  
  cat("\nCoefficients:\n")
  print(t1)
  l.RSS<-sum(gw.lm.out$LM_GofFit$LM_Res^2)
  mean.y<-mean(gw.lm.out$LM_GofFit$y)
  TSS<-sum((gw.lm.out$LM_GofFit$y-mean.y)^2)
  

  l.r<-1-(l.RSS/TSS)
  l.r.adj<-1-((l.RSS/(Obs-ModelVarNo-1))/(TSS/(Obs-1)))
  cat("\nResidual Sum of Squares:", l.RSS)
  cat("\nR-squared:", l.r)
  cat("\nAdjusted R-squared:", l.r.adj)
  
  return(gw.lm.out)
}