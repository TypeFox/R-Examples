
DensParcorr <- function(data,
                        select = FALSE,
                        dens.level = "plateau",
                        plateau.thresh = 0.01,
                        Parcorr.est = NULL,
                        directory = NULL
                        )
{ 
  if(is.null(directory)) directory = paste(getwd(),"/DensParcorr.output",sep="")
  
  lambda.max=0.6
  lambda.min=1e-8
  nlambda=10
  lambda = 10^(seq(log10(lambda.min), log10(lambda.max),length.out = nlambda))
  
  #### Calculate Precision Matrix ####
  data = as.matrix(data)

  if(!is.null(Parcorr.est))
  {
    if(is.null(Parcorr.est$precision.list)|is.null(Parcorr.est$lambda)) 
      stop("Parcorr.est is modified!")
  
    Prec.mat = Parcorr.est$precision.list
    lambda = Parcorr.est$lambda
    
    Prec.mat = Prec.mat[order(lambda)]
    lambda = lambda[order(lambda)]
    
  }else
  {
    lambda = lambda[order(lambda)]
    Prec.mat = clime(data,lambda = lambda)$Omegalist
  }

  # Calculate Dens Values for Precision Matrix
  dens = vector()
  for(i in 1:length(Prec.mat))
  {
    dens[i]=prec2dens(Prec.mat[[i]])
  }
      
  # To guarantee the maximum of Dens is valid
  while(abs(1 - sort(dens)[length(dens)-1]/max(dens))>0.05)
  {
    lambda = c(min(lambda)/10,lambda)
    Prec.mat = append(clime(data,lambda = lambda[1])$Omegalist,Prec.mat)
    dens = c(prec2dens(Prec.mat[[1]]),dens)
  }

  dir.create(directory)

  #### Based on different Tuning Parameter Selection Method ####
  if(select)
  {
      
      if(dens.level=="plateau")
      {
        select.index = max(which((1-dens/max(dens))<=plateau.thresh))
        
        while(dens[select.index]==max(dens))
        {
          lam.min = lambda[select.index]
          lam.max = lambda[select.index+1] 
          lambda2 = 10^(seq(log10(lam.min), log10(lam.max),length.out = 4))[2:3]
          Prec.mat2 = clime(data,lambda = lambda2)$Omegalist
          Prec.mat = append(Prec.mat2,Prec.mat)
          lambda = c(lambda2,lambda)
          
          Prec.mat = Prec.mat[order(lambda)]
          lambda = lambda[order(lambda)]
          for(i in 1:length(Prec.mat))
          {
            dens[i]=prec2dens(Prec.mat[[i]])
          }
          select.index = max(which((1-dens/max(dens))<=plateau.thresh))
        }

      }else if(is.numeric(dens.level)&dens.level<1&dens.level>0)
      {
        select.index = which.min(abs(dens-max(dens)*dens.level))
        
        while(abs(dens.level - dens[select.index]/max(dens))>0.05)
        {
          if(max(dens)*dens.level > dens[select.index])
          {
            lam.min = lambda[select.index-1]
            lam.max = lambda[select.index] 
          }else
          {
            lam.min = lambda[select.index]
            lam.max = lambda[select.index+1] 
          }
          lambda2 = 10^(seq(log10(lam.min), log10(lam.max),length.out = 4))[2:3]
          
          Prec.mat2 = clime(data,lambda = lambda2)$Omegalist
          Prec.mat = append(Prec.mat2,Prec.mat)
          lambda = c(lambda2,lambda)
          
          Prec.mat = Prec.mat[order(lambda)]
          lambda = lambda[order(lambda)]
          for(i in 1:length(Prec.mat))
          {
            dens[i]=prec2dens(Prec.mat[[i]])
          }
          select.index = which.min(abs(dens-max(dens)*dens.level))
        }
      }
      
      log10lam = -log(lambda)/log(10)
      
      
      
      png(filename = paste(directory,"/Dens.trace.plot.png",sep=""))
        par(mfrow=c(1,1))
        plot(log10lam[order(log10lam)],dens[order(log10lam)]/max(dens),xlab="-log10(lambda)",ylab="Percentage Dens",type="b")
        points(-log(lambda[select.index])/log(10),(dens/max(dens))[select.index],col="red",pch=19)
        legend("bottomright",legend="Selected",col="red",pch=19)
      dev.off()
    
  }else
  {
      png(filename = paste(directory,"/Dens.trace.plot.png",sep=""))
        par(mfrow=c(1,1))
        plot(-log(lambda[order(-log(lambda)/log(10))])/log(10),dens[order(-log(lambda)/log(10))]/max(dens),xlab="-log10(lambda)",ylab="Percentage Dens",type="b")
      dev.off()
 
  }
  
  #### Calculate Partial Correlation Matri from Precision Matrix ####
  
  dir.create(paste(directory,"/partial.correlation.matrix",sep=""))
  dir.create(paste(directory,"/precision.matrix",sep=""))
  
  Par.mat = list()
  for(i in 1:length(lambda))
  {
    Par.mat[[i]] = prec2part(Prec.mat[[i]])
    
       png(filename = paste(directory,"/partial.correlation.matrix/",i,".png",sep=""))
       heatmap.2(Par.mat[[i]],Rowv=F,Colv=F,scale="none",trace="none",density.info="none",xlab="",ylab="",
              main=paste("Dens Percentage=",round(dens[i]/max(dens),3),sep=""),col=bluered,dendrogram="none")
       dev.off()
       write.table(Par.mat[[i]],file = 
                     paste(directory,"/partial.correlation.matrix/",i,"_Partial.DensPercentage_",
                           round(dens[i]/max(dens),3),".txt",sep=""),
                   col.names = F,row.names = F)
       
       png(filename = paste(directory,"/precision.matrix/",i,".png",sep=""))
       heatmap.2(Par.mat[[i]],Rowv=F,Colv=F,scale="none",trace="none",density.info="none",xlab="",ylab="",
                 main=paste("Dens Percentage=",round(dens[i]/max(dens),3),sep=""),col=bluered,dendrogram="none")
       dev.off()
       write.table(Prec.mat[[i]],file = 
                     paste(directory,"/precision.matrix/",i,"_Precision.DensPercentage_",
                           round(dens[i]/max(dens),3),".txt",sep=""),
       col.names = F,row.names = F)
  }
  
  #### Summarize the results ####
  Results = list()
  if(select)
  {
      Results$selected.partial.corr = Par.mat[[select.index]]
      Results$selected.precision = Prec.mat[[select.index]]
      Results$selected.lambda = lambda[[select.index]]
      
      Results$partial.corr.list = Par.mat
      Results$precision.list = Prec.mat
      Results$lambda.list = lambda
      
  }else
  {
      Results$partial.corr.list = Par.mat
      Results$precision.list = Prec.mat
      Results$lambda.list = lambda
      Results$Dens = dens
      Results$Dens.Percentage = round(dens/max(dens),3)
  }
  
  if(select)
  {
      if(dens.level=="plateau")
      {
        Results$selection.method = "Dens-plateau"
        
      }else if(is.numeric(dens.level)&dens.level<1&dens.level>0)
      {
        Results$selection.method =  paste(round(dens.level*100,1),"% Dens (Actual=",round(dens[select.index]/max(dens)*100,1),"%)",sep="")
      }
      Results$Dens = dens
      Results$Dens.Percentage = round(dens/max(dens),3)
  }
  
  print("Figures, Estimated Precision Matrices and Partial Correlation are outputed in ")
  print(directory)
  
  return(Results)
}
