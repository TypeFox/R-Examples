gw.glm<-function(formula, family, dframe, bw, kernel, coords)
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
    
    Gl.Model <- eval(substitute(glm(formula, family=family, data = dframe)))
    
    RNames <- names(Gl.Model$coefficients)
    PValueNames <- paste("P",RNames, sep="_")

    ModelVarNo <- length(Gl.Model$coefficients)
    
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

    GGLM_LEst <- as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), RNames[1:ModelVarNo]))
    GGLM_LPvalues <- as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), PValueNames[1:ModelVarNo]))
    GGLM_GofFit <- data.frame(y=numeric(0), GLM_yfit=numeric(0), GLM_Res=numeric(0), 
                            GLM_AIC=numeric(0), GLM_Deviance=numeric(0),GLM_DevExpl=numeric(0),
                            GLM_Overdispersion=numeric(0),GLM_OverDis1=numeric(0))
  
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
    
    #Calculate WGLM
    Lcl.Model<-eval(substitute(glm(formula, family=family, data = SubSet, weights=Wts)))
    
    #Store in table
    GGLM_LEst[m,]<-Lcl.Model$coefficients
    GGLM_LPvalues[m,]<-summary(Lcl.Model)$coefficients[,4]
    
    GGLM_GofFit[m,1]<-Gl.Model$y[m]
    GGLM_GofFit[m,2]<-exp(sum(Gl.Model$model[m,2:ModelVarNo] * Lcl.Model$coefficients[2:ModelVarNo]) + sum(Lcl.Model$coefficients[1]))
   
    GGLM_GofFit[m,4]<-AIC(Lcl.Model)
    GGLM_GofFit[m,5]<-deviance(Lcl.Model)
    GGLM_GofFit[m,6]<-100*((Lcl.Model$null.deviance-Lcl.Model$deviance)/Lcl.Model$null.deviance)
    GGLM_GofFit[m,7]<-Lcl.Model$deviance/Lcl.Model$df.residual
    GGLM_GofFit[m,8]<-var(Lcl.Model$model[,1])/mean(Lcl.Model$model[,1])
  
  }
   
    Adj.f.l <- sum(GGLM_GofFit[,1])/sum(GGLM_GofFit[,2])
    GGLM_GofFit[,3] <- GGLM_GofFit[,1] - Adj.f.l*GGLM_GofFit[,2]
    
  gw.glm.out<-list(GGLM_LEst=GGLM_LEst, GGLM_LPvalues=GGLM_LPvalues, GGLM_GofFit=GGLM_GofFit)
  
  cat("\n--------------- Local Model Summary ---------------\n")
  
  cat("\nResiduals:\n")
  print(summary(gw.glm.out$GGLM_GofFit$GLM_Res))
  
  t1 <- data.frame(Min = apply(gw.glm.out$GGLM_LEst, 2, min), Max = apply(gw.glm.out$GGLM_LEst, 2, max), 
            Mean = apply(gw.glm.out$GGLM_LEst, 2, mean), StD = apply(gw.glm.out$GGLM_LEst, 2, sd))
  
  cat("\nCoefficients:\n")
  print(t1)

  l.RSS<-sum(gw.glm.out$GGLM_GofFit$GLM_Res^2)
  
  l.r2<-1-(l.RSS/TSS)
  l.r2.adj<-1-((l.RSS/(Obs-ModelVarNo-1))/(TSS/(Obs-1)))
  cat("\nResidual Sum of Squares:", l.RSS)
  cat("\nR-squared:", l.r2)
  cat("\nAdjusted R-squared:", l.r2.adj)
  
  return(gw.glm.out)
}