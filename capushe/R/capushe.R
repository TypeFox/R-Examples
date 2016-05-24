setClass(
  Class="Capushe",
  representation=representation(
    DDSE="DDSE",
    Djump="Djump",
    AIC_capushe="list",
    BIC_capushe="list",
    n="numeric")
)


setMethod("print","Capushe",
          function(x,...){
            result=list(x@DDSE@model,x@Djump@model)
            names(result)=c("DDSE_model","Djump_model")
            if (x@n!=0){
              Tempo=list(x@AIC_capushe$model,x@BIC_capushe$model)
              names(Tempo)=c("AIC_model","BIC_model")
              result=c(result,Tempo)
            }
            print(result)
          }
)

setMethod("show","Capushe",
          function(object){
            result=list(object@DDSE@model,object@Djump@model)
            names(result)=c("DDSE_model","Djump_model")
            if (object@n!=0){
              Tempo=list(object@AIC_capushe$model,object@BIC_capushe$model)
              names(Tempo)=c("AIC_model","BIC_model")
              result=c(result,Tempo)
            }
            print(result)
          }
)





capushe = function(data,n=0,pct=0.15,point=0,psi.rlm=psi.bisquare,scoef=2,Careajump=0,Ctresh=0){
  
  if(!(is.numeric(n))){
    stop("n must be numeric")
  }
  if(!(is.numeric(pct))){
    stop("pct must be numeric")
  }
  if(!(is.numeric(point))){
    stop("point must be numeric")
  }
  if (!(floor(point)==point)|(point<0)){
    stop("point must be an positive integer")
  }
  if(!(is.numeric(scoef))){
    stop("scoef must be numeric")
  }
  
  if (is.character(data)){data=read.table(data)}
  
  if (length(data[1,])!=4){
    stop("data must have 4 columns with name, penshape, complexity and contrast")
  }
  data=data.frame(data)
  names(data)=c("model","pen","complexity","contrast")
  
  if(any(is.na(data))){
    bad=which(is.na(data),arr.ind=TRUE)[,1]
    repbad=which(duplicated(bad),arr.ind=TRUE)
    if (length(repbad)!=0){bad=bad[-repbad]}
    data=data[-bad,]
    if (length(bad)==1){
      warning(paste("1 line has been removed by capushe"))
    }
    else{warning(paste(as.character(length(bad)),"lines have been removed by capushe"))}
  }
  mlength=length(data$model)
  if (mlength<10){
    stop("At least 10 observations are needed")
  }
  
  if (!is.numeric(data[,2])){
    stop("pen must be numeric")
  }
  if (!is.numeric(data[,3])){
    stop("complexity must be numeric")
  }
  if (!is.numeric(data[,4])){
    stop("contrast must be numeric")
  }
  
  if ((length(data$pen)!=mlength)|(length(data$complexity)!=mlength)|(length(data$contrast)!=mlength)){
    stop("The lengths of the columns are not equal")
  }
  if ((pct<0)|(pct>1)){
    stop("pct must be between 0 and 1")
  }
  if (point>mlength){
    stop("point must be smaller than the number of models")
  }
  if (!prod(data$complexity>=0)){
    stop("Complexity must be positive")
  }
  if (!(floor(n)==n)|(n<0)){
    stop("The number of observations must be an positive integer")
  }
  if(n==Inf){
    stop("n must be different of Inf")
  }
  if (scoef<1){
    stop("scoef must be greater than one")
  }
  if(scoef==Inf){
    stop("scoef must be different of Inf")
  }
  if(!(is.numeric(Careajump))){
    stop("Caerajump must be numeric")
  }
  if(Careajump<0){
    stop("Careajump must be positive")
  }
  if(Careajump>=1){
    stop("Careajump must be less than one")
  }
  if(!(is.numeric(Ctresh))){
    stop("Ctresh must be numeric")
  }
  if(Ctresh<0){
    stop("Ctresh must be positive") 
  }
  if(Ctresh>max(data$complexity)){
    stop("Ctresh must be less than the maximum complexity")
  }
  
  
  slope=DDSE(data,pct=pct,point=point,psi.rlm=psi.rlm,scoef=scoef)
  if ((Careajump==0)*(n!=0)){
    jump=Djump(data,scoef=scoef,Careajump=sqrt(log(n)/n),Ctresh=Ctresh)
  }else{
    jump=Djump(data,scoef=scoef,Careajump=Careajump,Ctresh=Ctresh)
  }
  
  if (slope@model!=jump@model){
    warning("Models selected by DDSE and Djump are different")
  }
  
  if (n!=0){
    aicapushe=AICcapushe(data,n)
    bicapushe=BICcapushe(data,n)
    new(Class="Capushe",DDSE=slope,Djump=jump,AIC_capushe=aicapushe,BIC_capushe=bicapushe,n=n)
  }
  else{
    new(Class="Capushe",DDSE=slope,Djump=jump,n=0)
  }
}


setMethod(
  f="validation",
  signature="Capushe",
  definition=function(x,data2,newwindow=TRUE,...){
    validation(x@DDSE,data2,newwindow)
  }
)

setMethod(
  f="plot",
  signature="Capushe",
  definition=function(x,y,newwindow=TRUE,ask=TRUE,...){
    
    if(!is.logical(newwindow)){
      stop("newwindow must be logical")
    }
    if(!is.logical(ask)){
      stop("ask must be logical")
    }
    refgraph=list(ask=par()$ask,
                  mgp=par()$mgp,
                  oma=par()$oma,
                  xaxs=par()$xaxs,
                  mfrow=par()$mfrow,
                  cex.axis=par()$cex.axis)
    
    plength=length(x@DDSE@graph$model)-1
    p=plength-x@DDSE@interval$point_using+2
    scoef=x@Djump@ModelHat$Kopt/x@Djump@ModelHat$kappa[x@Djump@ModelHat$JumpMax+1]
    
    Model=x@DDSE@ModelHat$model_hat[x@DDSE@ModelHat$point_breaking[x@DDSE@ModelHat$imax]]
    
    if (newwindow){	dev.new(width=14)}
    for (i in 1:(newwindow+1)){
      par(ask=ask,mgp = c(3, 0.5, 0))
      layout(matrix(c(1,1,1,1,1,1,1,1,2,3,2,3,2,3),ncol=7))
      Couleur=c(rep("blue",(p-1)),rep("red",(plength-p+2)))
      plot(x=x@DDSE@graph$pen,y=-x@DDSE@graph$contrast,main="Contrast representation",xlab=expression(paste(pen[shape](m)," (labels : Model names)")),ylab=expression(-gamma[n] (hat(s)[m])),col=Couleur,xaxt = "n")
      
      legend("bottomright",legend=c(paste("The regression line is computed with",as.character(plength-p+2),"points")),lty=1,col=c("red"))
      
      labelpen=x@DDSE@graph$pen
      labelmodel=x@DDSE@graph$model
      pas=(max(labelpen)-min(labelpen))/50*11.5/par("din")[1]
      leng=length(labelpen)
      Coord=c()
      i=1
      j=2
      while (j<leng){
        while ((j<leng)*((labelpen[j]-labelpen[i])<pas)){
          Coord=c(Coord,-j)
          j=j+1
        }
        i=j
        j=j+1
      }
      if (length(Coord)>0){
        labelpen=labelpen[Coord]
        labelmodel=labelmodel[Coord]
      }
      axis(1,labelpen,labelmodel,las=2)
      
      Cof=x@DDSE@graph$reg$coefficients
      abline(a=Cof[1],b=Cof[2],col="red")
      abscisse=seq(plength+1,2)
      plot(x@DDSE@kappa,main="Successive slope values",xlab=expression(paste("Number of points (",pen[shape](m),",",-gamma[n] (hat(s)[m]),") for the regression")),ylab=expression(kappa),type="l",xaxt = "n")
      
      abscisse2=seq(1,plength+1,length=floor(13*11.5/par("din")[1]))
      axis(1,abscisse2,abscisse[abscisse2])
      
      points(x@DDSE@kappa,col="red",pch=18)
      plot(x@DDSE@ModelHat$model_hat,main="Selected models with respect to the successive slope values",xlab=expression(paste("Number of points (",pen[shape](m),",",-gamma[n] (hat(s)[n]),") for the regression")),ylab="Model",pch=20,yaxt="n",xaxt = "n")
      
      axis(1,abscisse2,abscisse[abscisse2])
      
      labelhat=x@DDSE@ModelHat$model_hat
      pas=(max(labelhat)-min(labelhat))/61*5.75/par("din")[2]
      leng=length(labelhat)
      Coord=c()
      i=1
      j=2
      while (j<leng){
        while ((j<leng)*((labelhat[j]-labelhat[i])<pas)){
          Coord=c(Coord,-j)
          j=j+1		
        }
        i=j
        j=j+1
      }
      Coord=c(Coord,-which(duplicated(labelhat)==TRUE))
      if (length(Coord)>0){
        labelhat=labelhat[Coord]
      }
      
      axis(2,labelhat,x@DDSE@graph$model[labelhat],las=1)
      
      lines(x=x@DDSE@ModelHat$point_breaking[x@DDSE@ModelHat$imax]:(x@DDSE@ModelHat$point_breaking[x@DDSE@ModelHat$imax]+x@DDSE@ModelHat$number_plateau[x@DDSE@ModelHat$imax]-1),y=rep(Model,x@DDSE@ModelHat$number_plateau[x@DDSE@ModelHat$imax]),col="blue")
      
    }
    
    
    ### Plot Djump
    leng=length(x@Djump@graph$model)
    mleng=length(x@Djump@ModelHat$model_hat)
    
    
    if (newwindow){	dev.new(width=15)}
    else{layout(c(1))}
    
    par(oma=c(0,3,0,0),xaxs="i")
    
    Intervalslope=scoef*x@DDSE@interval$interval
    if (length(x@Djump@Area)==4){
      Absmax=max((scoef+1)/scoef*x@Djump@ModelHat$Kopt,5/4*Intervalslope[2],5/4*scoef*x@Djump@Area$kappasup)
    }else{
      Absmax=max((scoef+1)/scoef*x@Djump@ModelHat$Kopt,5/4*Intervalslope[2])
    }
    
    Ordmin=max(which((x@Djump@ModelHat$kappa<=Absmax)==TRUE))
    plot(x=x@Djump@ModelHat$Kopt,y=x@Djump@graph$complexity[x@Djump@graph$Modopt],xlim=c(0,Absmax),ylim=c(x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[Ordmin]],x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[1]]),ylab="",xlab=expression(paste("Values of the penalty constant ",kappa)),yaxt="n",pch=4,col="blue",main="Dimension Jump",lwd=3,xaxt="n")
    complex=x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[1:Ordmin]]
    complexgreen=complex
    for (i in 1:(mleng-1)){
      lines(x=c(x@Djump@ModelHat$kappa[i],x@Djump@ModelHat$kappa[i+1]),y=c(complex[i],complex[i]))
    }
    if (x@Djump@ModelHat$kappa[i+1]<Absmax){
      lines(x=c(x@Djump@ModelHat$kappa[mleng],Absmax),y=c(complex[mleng],complex[mleng]))
    }
    lines(x=c(x@Djump@ModelHat$Kopt/scoef,x@Djump@ModelHat$Kopt/scoef),y=c(complex[x@Djump@ModelHat$JumpMax],complex[x@Djump@ModelHat$JumpMax+1]),col="blue",lty=2)
    ordon=paste(as.character(complex),"(",as.character(x@Djump@graph$model[x@Djump@ModelHat$model_hat[1:Ordmin]]),")")
    par(cex.axis=0.6)
    
    complex2=complex
    pas=(max(complex2)-min(complex2))/39*5.6/par("din")[2]
    leng2=length(complex2)
    Coord=c()
    i=leng2
    j=leng2-1
    while (j>0){
      while ((j>1)*((complex2[j]-complex2[i])<pas)){
        Coord=c(Coord,-j)
        j=j-1		
      }
      if ((j==1)*((complex2[j]-complex2[i])<pas)){
        Coord=c(Coord,-j)
        j=j-1		
      }
      i=j
      j=j-1
    }
    if (x@Djump@graph$Ctresh>0){
      Coord=c(Coord,-which(abs(complex-x@Djump@graph$Ctresh)<pas))
      lines(c(par("usr")[1],par("usr")[2]),c(x@Djump@graph$Ctresh,x@Djump@graph$Ctresh),lty=4,col="chocolate1")
      axis(2,x@Djump@graph$Ctresh,expression(C[tresh]),las=2)
    }
    if (length(Coord)>0){
      complex2=complex2[Coord]
      ordon=ordon[Coord]
    }
    axis(2,complex2,ordon,las=2)
    
    
    par(cex.axis=1)
    axis(2,labels="Selected Complexity (Model)",outer=TRUE,at=(x@Djump@graph$complexity[leng]+x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[Ordmin]])/2,lty=0)
    if (scoef==1){
      axis(1,labels=  expression(kappa[opt] == hat(kappa)^{dj}),at=x@Djump@ModelHat$Kopt)
    }
    else {  
      axis(1,labels=expression(hat(kappa)^{dj}),at=x@Djump@ModelHat$Kopt/scoef)
      axis(1,labels=  expression(kappa[opt]),at=x@Djump@ModelHat$Kopt)
    }
    
    if (scoef<=2){
      absci=seq(0,Absmax,by=x@Djump@ModelHat$Kopt/scoef/3)[c(-4)]
      Empty=c(which(abs(absci-x@Djump@ModelHat$Kopt)<absci[2]/4*(12.093749/par("din")[1])^2*scoef),which(abs(absci-x@Djump@ModelHat$Kopt/scoef)<absci[2]/4*(12.093749/par("din")[1])^2*scoef))
    }
    if ((2<scoef)*(scoef<=4)){
      absci=seq(0,Absmax,by=x@Djump@ModelHat$Kopt/scoef/2)[c(-3)]
      Empty=c(which(abs(absci-x@Djump@ModelHat$Kopt)<absci[2]/15*(12.093749/par("din")[1])^2*scoef),which(abs(absci-x@Djump@ModelHat$Kopt/scoef)<absci[2]/15*(12.093749/par("din")[1])^2*scoef))
    }
    if (scoef>4){
      absci=seq(0,Absmax,by=x@Djump@ModelHat$Kopt/scoef)[c(-2)]
      Empty=c(which(abs(absci-x@Djump@ModelHat$Kopt)<absci[2]/30*(12.093749/par("din")[1])^2*scoef),which(abs(absci-x@Djump@ModelHat$Kopt/scoef)<absci[2]/30*(12.093749/par("din")[1])^2*scoef))
    }
    
    if (length(Empty)>0){
      absci=absci[-Empty]
    }
    axis(1,absci,signif(absci,2))
    
    redmin=max(which(x@Djump@ModelHat$kappa<Intervalslope[1]))
    redmax=max(which(x@Djump@ModelHat$kappa<Intervalslope[2]))
    
    # 	if (redmin==(redmax-1)){
    # 	  lines(x=c(Intervalslope[1],Intervalslope[2]),y=c(complex[redmin],complex[redmin]),col="red",lwd=2)
    # 	}else{
    # 	  lines(x=c(Intervalslope[1],x@Djump@ModelHat$kappa[redmin+1]),y=c(complex[redmin],complex[redmin]),col="red",lwd=2)
    # 	  if ((redmin+1)<(redmax-1)){
    # 	    for (i in (redmin+1):(redmax-2)){
    # 	      lines(x=c(x@Djump@ModelHat$kappa[i],x@Djump@ModelHat$kappa[i+1]),y=c(complex[i],complex[i]),col="red",lwd=2)
    # 	    }
    # 	  }
    # 	  lines(x=c(x@Djump@ModelHat$kappa[redmax-1],Intervalslope[2]),y=c(complex[redmax-1],complex[redmax-1]),col="red",lwd=2)
    # 	}
    if (redmin==redmax){
      lines(x=c(Intervalslope[1],Intervalslope[2]),y=c(complex[redmin],complex[redmin]),col="red",lwd=2)
    }else{
      R=redmin:redmax
      Kap=c(Intervalslope[1],x@Djump@ModelHat$kappa[R[-1]],Intervalslope[2])
      for (r in 1:length(R)){
        lines(x=Kap[c(r,r+1)],y=rep(complex[R[r]],2),col="red",lwd=2)
      }
    }
    
    if (length(x@Djump@Area)==4){
      lines(c(x@Djump@Area$kappainf,x@Djump@Area$kappasup),c(par("usr")[3],par("usr")[3]),col="green4",lwd=2)
      yarea=(x@Djump@graph$complexity[x@Djump@ModelHat$model_hat[Ordmin]]+par("usr")[3])/2
      lines(c(x@Djump@Area$kappainf,x@Djump@Area$kappainf),c(par("usr")[3],yarea),col="green4",lwd=2)
      lines(c(x@Djump@Area$kappasup,x@Djump@Area$kappasup),c(par("usr")[3],yarea),col="green4",lwd=2)
      
      greenmin=max(which(x@Djump@ModelHat$kappa<scoef*x@Djump@Area$kappainf))
      greenmax=max(which(x@Djump@ModelHat$kappa<scoef*x@Djump@Area$kappasup))
      
#       G=unique(c(scoef*x@Djump@Area$kappainf,
#           x@Djump@ModelHat$kappa[(scoef*x@Djump@Area$kappainf<x@Djump@ModelHat$kappa)&&(scoef*x@Djump@Area$kappasup>x@Djump@ModelHat$kappa)],
#           scoef*x@Djump@Area$kappasup))
      
      if (greenmin==greenmax){
        lines(x=scoef*c(x@Djump@Area$kappainf,x@Djump@Area$kappasup),y=c(complexgreen[greenmin],complexgreen[greenmin]),col="green4",lwd=4,lty=3)
      }else{
        G=greenmin:greenmax
        Kap=c(scoef*x@Djump@Area$kappainf,x@Djump@ModelHat$kappa[G[-1]],scoef*x@Djump@Area$kappasup)
        for (g in 1:length(G)){
          lines(x=Kap[c(g,g+1)],y=rep(complexgreen[G[g]],2),col="green4",lwd=4,lty=3)
        }
      }
      
      legend("topright",legend=c("Maximal Jump",paste(as.character(scoef), "x Slope estimation interval"),"Jump area",paste(as.character(scoef), "x Jump area")),lty=c(2,1,1,3),col=c("blue","red","green4","green4"),lwd=c(1,2,2,4),bg = "white")
    }else{
      legend("topright",legend=c("Maximal Jump",paste(as.character(scoef), "x Slope estimation interval")),lty=c(2,1),col=c("blue","red"),lwd=c(1,2),bg = "white")
    }
    
    par(refgraph)
  }
)


setMethod("summary","Capushe",
          function(object){
            
            cat("\nModel selected by DDSE: ",as.character(object@DDSE@model),"\n",sep="")
            
            cat("\nCorresponding Slope interval: [",as.character(signif(object@DDSE@interval$interval[1],3)),";",as.character(signif(object@DDSE@interval$interval[2],3)),"](",as.character(signif(object@DDSE@interval$percent_of_points,3)*100),"% of points)\n",sep="")
            cat("Selected model complexity: ",as.character(object@DDSE@graph$Complexity),"\n\n",sep="")
            
            
            cat("\nModel selected by Djump: ",as.character(object@Djump@model),"\n",sep="")
            
            cat("\nEstimated penalty constant: ",as.character(signif(object@Djump@ModelHat$kappa[object@Djump@ModelHat$JumpMax+1],3)),"\n",sep="")
            cat("Selected model complexity: ",as.character(object@Djump@graph$complexity[object@Djump@graph$Modopt]),"\n\n",sep="")
            if (length(object@Djump@Area)==4){
              cat("Jump area: [",as.character(signif(object@Djump@Area$kappainf,3)),";",as.character(signif(object@Djump@Area$kappasup,3)),"]\n",sep="")
              if (length(object@Djump@Area$Modeloptarea)>1){
                cat("Models selected by jump area: ",sep="")
                cat(cat(as.character(object@Djump@graph$model[sort(object@Djump@Area$Modeloptarea)]),sep=", "),"\n\n",sep="")
              }else{
                cat("Model selected by jump area: ",as.character(object@Djump@graph$model[object@Djump@Area$Modeloptarea]),"\n\n",sep="")
              }
            }
            
            if (object@n!=0){
              cat("\nModel selected by AIC: ",as.character(object@AIC_capushe$model),"\n",sep="")
              cat("\nModel selected by BIC: ",as.character(object@BIC_capushe$model),"\n",sep="")
            }
          }
)