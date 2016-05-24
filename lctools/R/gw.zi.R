gw.zi <- function(formula, family, dframe, bw, kernel, coords)
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

  Gl.Model<-eval(substitute(zeroinfl(formula, data = dframe, dist=family)))
  
  RNames<-names(Gl.Model$coefficients$count)
  CountNames<-paste("CM",RNames, sep="_")
  ZeroNames<-paste("ZM",RNames, sep="_")
  PValueNames<-paste("P",RNames, sep="_")
  PValueNamesCount<-paste("P_CM",RNames, sep="_")
  PValueNamesZero<-paste("P_ZM",RNames, sep="_")
  
  ModelVarNo<-length(Gl.Model$coefficients$count)
  
  cat("\nNumber of Variables:", ModelVarNo-1)
  cat("\n--------------- Global Model Summary ---------------\n")
  print(summary(Gl.Model))
  
  Gl.Model.res<-Gl.Model$y-Gl.Model$fitted.values
  Gl.RSS<-sum(Gl.Model.res^2)
  mean.y<-mean(Gl.Model$y)
  TSS<-sum((Gl.Model$y-mean.y)^2)
  
  Gl.r2<-1-(Gl.RSS/TSS)
  Gl.r2.adj<-1-((Gl.RSS/(Obs-ModelVarNo-1))/(TSS/(Obs-1)))
  cat("\nResidual Sum of Squares:", Gl.RSS)
  cat("\nR-squared:", Gl.r2)
  cat("\nAdjusted R-squared:", Gl.r2.adj)
  
  DistanceT <- dist(coords)
  Dij <- as.matrix(DistanceT)
  
  ZI_LEst_count<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), CountNames[1:ModelVarNo]))
  ZI_LEst_zero<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), ZeroNames[1:ModelVarNo]))
  ZI_LPvalues_count<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), PValueNamesCount[1:ModelVarNo]))
  ZI_LPvalues_zero<-as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), PValueNamesZero[1:ModelVarNo]))
  ZI_GofFit<-data.frame(y=numeric(0), ZI_yfit=numeric(0), ZI_Res=numeric(0), ZI_AIC=numeric(0), ZI_LogLik=numeric(0),ZI_OverDis1=numeric(0))
  
  for(m in 1:Obs){

    #Get the data
    DNeighbour <- Dij[,m]
    DataSet <- data.frame(dframe, DNeighbour=DNeighbour)

    #Sort by distance
    DataSetSorted <- DataSet[order(DataSet$DNeighbour),]
    
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
    ZI_LPvalues_count[m,]<-summary(ZI.Model)$coefficients$count[,4]  
    ZI_LPvalues_zero[m,]<-summary(ZI.Model)$coefficients$zero[,4]  
    
    ZI_GofFit[m,1]<-Gl.Model$y[m]
    ZI_GofFit[m,2]<-exp(sum(Gl.Model$model[m,2:ModelVarNo] * ZI.Model$coefficients$count[2:ModelVarNo]) + sum(ZI.Model$coefficients$count[1]))
    # ZI_GofFit[m,3]<-ZI_GofFit[m,1]-ZI_GofFit[m,2]
    
    ZI_GofFit[m,4]<-AIC(ZI.Model)
    ZI_GofFit[m,5]<-ZI.Model$loglik
    
    ZI_GofFit[m,6]<-var(ZI.Model$model[,1])/mean(ZI.Model$model[,1])
    
    #cat("-",m)
  }
  
  Adj.f.l <- sum(ZI_GofFit[,1])/sum(ZI_GofFit[,2])
  ZI_GofFit[,3] <- ZI_GofFit[,1] - Adj.f.l*ZI_GofFit[,2]
  
  gw.zi.results<-list(ZI_LEst_count=ZI_LEst_count, ZI_LEst_zero=ZI_LEst_zero, ZI_LPvalues_count=ZI_LPvalues_count, ZI_LPvalues_zero=ZI_LPvalues_zero, ZI_GofFit=ZI_GofFit)

  cat("\n--------------- Local Model Summary ---------------\n")
  
  cat("\nResiduals:\n")
  print(summary(gw.zi.results$ZI_GofFit$ZI_Res))
  
  t1 <- data.frame(Min = apply(gw.zi.results$ZI_LEst_count, 2, min), Max = apply(gw.zi.results$ZI_LEst_count, 2, max), 
                   Mean = apply(gw.zi.results$ZI_LEst_count, 2, mean), StD = apply(gw.zi.results$ZI_LEst_count, 2, sd))
  t2 <- data.frame(Min = apply(gw.zi.results$ZI_LEst_zero, 2, min), Max = apply(gw.zi.results$ZI_LEst_zero, 2, max), 
                   Mean = apply(gw.zi.results$ZI_LEst_zero, 2, mean), StD = apply(gw.zi.results$ZI_LEst_zero, 2, sd))
  
  cat("\nCoefficients:\n")
  cat("\nCount Model:\n")
  print(t1)
  cat("\nZero Model:\n")
  print(t2)
  
   
  l.RSS<-sum(gw.zi.results$ZI_GofFit$ZI_Res^2)
  
  l.r2<-1-(l.RSS/TSS)
  l.r2.adj<-1-((l.RSS/(Obs-ModelVarNo-1))/(TSS/(Obs-1)))
  cat("\nResidual Sum of Squares:", l.RSS)
  cat("\nR-squared:", l.r2)
  cat("\nAdjusted R-squared:", l.r2.adj)
  
  return(gw.zi.results)
}