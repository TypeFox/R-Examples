plotMeanSD <- function(x,y,y.se,ylabel,xlabel){
      ly <- y -y.se
      uy <- y+y.se
      plot(x,y,xlim=range(x),ylim=range(c(y,ly,uy)),xlab=xlabel,ylab=ylabel,cex.axis=1.5,cex=1.5,cex.lab=1.5,type="n",xaxt="n")
      points(x,y,pch=19,cex.axis=1.5,cex=1.5,cex.lab=1.5)
      axis(1, at = c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95), labels = c(0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95),cex.axis=1.5,cex=1.5,cex.lab=1.5)
      for(i in 1:length(x)){
         lines(c(x[i],x[i]),c(ly[i],uy[i]),cex.axis=1.5,cex=1.5,cex.lab=1.5)
         lines(c(x[i]-0.01,x[i]+0.01),c(ly[i],ly[i]),cex.axis=1.5,cex=1.5,cex.lab=1.5)
         lines(c(x[i]-0.01,x[i]+0.01),c(uy[i],uy[i]),cex.axis=1.5,cex=1.5,cex.lab=1.5)

      }

}


plotLambda <- function(lambdaChoiceOutput,output=c("wss","ncluster","pwss","ch","h")){
  nGene=100
  if(length(output)>1){
  warning("You can Specify only one value at a time. The output of the first value will be provided ")
  output <- output[1]
  }
  Lambda <- lambdaChoiceOutput[[1]][,1]
  selectOutput <- c("wss","ncluster","pwss","ch","h")%in% output
  OutputIndex <- c(1:5)[selectOutput]
  if(OutputIndex==1){ ## WSS
        wOutput <- NULL
       wssMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  wOutput <- cbind(wOutput,lambdaChoiceOutput[[x]][,2]))
       wssMean <- rowMeans(wssMatrix)
       wssSe <- apply(wssMatrix,1,function(x) sd(x)/sqrt(length(x)))
       plotMeanSD(x=Lambda,y=wssMean,y.se=wssSe,ylabel="WSS",xlabel="Lambda")
  }

 if(OutputIndex==2){ ##ncluster
       nOutput <- NULL
       ncMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  nOutput <- cbind(nOutput,lambdaChoiceOutput[[x]][,4]))
       ncMean <- rowMeans(ncMatrix)
       ncSe <- apply(ncMatrix,1,function(x) sd(x)/sqrt(length(x)))
       plotMeanSD(x=Lambda,y=ncMean,y.se=ncSe,ylabel="No. of Clusters",xlabel="Lambda")
  }

 if(OutputIndex==3){ ##pwss
       nOutput <- wOutput <- NULL
       ncMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  nOutput <- cbind(nOutput,lambdaChoiceOutput[[x]][,4]))
       wssMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  wOutput <- cbind(wOutput,lambdaChoiceOutput[[x]][,2]))
       pwss <-  wssMatrix + 2*ncMatrix
       pwssMean <- rowMeans(pwss)
       pwssSe <- apply(pwss,1,function(x) sd(x)/sqrt(length(x)))
       plotMeanSD(x=Lambda,y=pwssMean,y.se=pwssSe,ylabel="pWSS",xlabel="Lambda")
  }


 if(OutputIndex==4){ ##CH
       nOutput <- wOutput <-tOutput <- NULL
       ncMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  nOutput <- cbind(nOutput,lambdaChoiceOutput[[x]][,4]))
       tssMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  tOutput <- cbind(tOutput,lambdaChoiceOutput[[x]][,3]))
       wssMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  wOutput <- cbind(wOutput,lambdaChoiceOutput[[x]][,2]))
       bssmatrix <- tssMatrix -  wssMatrix
       num <-  bssmatrix/(ncMatrix)
       den <- wssMatrix/(nGene - ncMatrix)
       ch <-  num/den
       chMean <- rowMeans(ch)
       chSe <- apply(ch,1,function(x) sd(x)/sqrt(length(x)))
       plotMeanSD(x=Lambda,y=chMean,y.se=chSe,ylabel="CH Values",xlabel="Lambda")
  }

 if(OutputIndex==5){ ##H
       nOutput <- wOutput <- NULL
       ncMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  nOutput <- cbind(nOutput,lambdaChoiceOutput[[x]][,4]))
       wssMatrix <- sapply(c(1:length(lambdaChoiceOutput)),function(x)  wOutput <- cbind(wOutput,lambdaChoiceOutput[[x]][,2]))
       wssMatrix1 <- wssMatrix[1:(length(Lambda)-1),]
       wssMatrix2 <- wssMatrix[2:length(Lambda),]
       ncMatrix2 <- ncMatrix[2:length(Lambda),]
       houtput <- ((wssMatrix1/wssMatrix2)-1)/(ncMatrix2)
       hMean <- rowMeans(houtput)
       hSe <- apply(houtput,1,function(x) sd(x)/sqrt(length(x)))
       plotMeanSD(x=Lambda[-length(Lambda)],y=hMean,y.se=hSe,ylabel="H Values",xlabel="Lambda")
  }

}


