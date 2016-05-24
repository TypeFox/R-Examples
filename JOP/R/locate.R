    ########################################################
    ####### Function that counts zeros after comma   #######
    ########################################################
    countdig<-function(x)
    {
      count<-0
      if(abs(x)<1e-8){
      return(0)}
      while(abs(x)<1)
      {
        x<-x*10
        count<-count+1
      }
      return(count)
    }
    ########################################################
    ####### Function checks if two numbers are equal #######
    ########################################################
    numeq<-function(x,y)
    {
      count<-0
      if(abs(x)<1e-8){
      return(abs(x-y)<1e-3)}
      x1<-x
      while(abs(x1)<1)
      {
        x1<-x1*10
        count<-count+1
      }
      tol<-1e-3*as.double(paste("1e-",count,sep=""))
      return(abs(x-y)<tol)
    }


#######################
#######################

locate <-
function(x,ncom=1,xlu=NULL,no.col=FALSE,standard=TRUE,col=1,lty=1,bty="l",las=1,adj=0.5,cex=1,cex.lab=1,cex.axis=1,xlab=c("Stretch Vector","Stretch Vector"),ylab=c("Parameter Setting","Predicted Response"),lwd=1,...)
{  
    ## Setting Values
    out<-x
    Wstart<-out$ValW[1]
    Wend<-out$ValW[2]  
    numbW<-dim(out$Responses)[1]
    nx<-dim(out$Parameters)[2]
    ny<-dim(out$Responses)[2]
    numbW<-dim(out$Responses)[1]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]  
    tau1<-out[[6]] 

    if(no.col==FALSE)
    {
      ## Setting Values
      cols<-1:(nx+ny)
      if((nx+ny)>=7)
      {
        cols<-1:(nx+ny+1)
        cols<-cols[-7]
      }
      ltys<-c(1:nx,1:ny)
    }
    if(no.col==TRUE)
    {
      col1<-gray(seq(0.2,0.5,length=nx))
      col2<-gray(seq(0.2,0.5,length=ny))
      cols<-c(col1,col2)
      ltys<-c(1:nx,1:ny)
    }
    if(no.col!=TRUE && no.col!=FALSE)
    {
      return(cat("no.col has to be logical!\n"))
    }
    if(length(col)==1)
    col<-cols
    if(length(lty)==1)
    lty<-ltys
    if(length(col)!=(nx+ny))
    col<-cols
    if(length(lty)!=(nx+ny))
    lty<-ltys
  
  
    if(!is.null(xlu))
    {
      if(any(xlu<1) || any(xlu>numbW))
      {
        cat("Choose a x-coordinate between 1 and numbW")
        cat("\n")
        return("Call locate again!")
      }
    }
    if(is.null(xlu))
    {
      cat("Choose your preferred point on the right plot please!\n")
      cat("\n")
      flush.console()
    }
  
    plot(x=out,no.col=no.col,standard=standard,col=col,lty=lty,bty=bty,las=las,adj=adj,cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,xlab=xlab,ylab=ylab,lwd=lwd,...) 
  

    ####
    ####  Location of points
    ####


    if(is.null(xlu))
    {
      if(!is.wholenumber(ncom))
      {
        return("Number of preferenced compromises has to be an integer!\n")
      }
      locp<-locator(n=ncom)
      xloc1<-round(locp$x)  
    }
    if(!is.null(xlu))
    {
      if(any(!is.wholenumber(xlu)))
      {
        return("The entries of xlu have to be integers!\n")
      }
      xloc1<-xlu  
    }
      
    ## Chosen Points:
    xlp<-xloc1
    optp<-optmatrix[xloc1,]
    reoptp<-reoptmatrix[xloc1,]
    ## Values to label the axes
    xaxis1<-1:numbW
    yaxis1<-seq(-round(max(abs(optmatrix)),digits=2+countdig(max(abs(optmatrix)))),round(max(abs(optmatrix)),digits=2+countdig(max(abs(optmatrix)))),round(max(abs(optmatrix)),digits=2+countdig(max(abs(optmatrix)))))

    optplot<-matrix(NaN,ncol=nx,nrow=numbW)
    for(i in 1:nx)
    {
      optplot[,i]<-optmatrix[,i]-min(optmatrix[,i])
      for(j in 1:numbW)
      { 
        optplot[j,i]<-ifelse(numeq(max(optmatrix[,i]),min(optmatrix[,i]))==TRUE,0.5,optplot[j,i]/(max(optmatrix[,i])-min(optmatrix[,i])))
      }  
    }
    
    par(mfrow=c(1,2),las=1,bty="l",adj=adj,mar=c(4,5,3,2))
    matplot(xaxis1,optplot,type="l",col=col[1:(nx)],lty=lty[1:(nx)],xlab=xlab[1],ylab="",ylim=c(0,1.25),xaxt="n",yaxt="n",cex.lab=cex.lab,lwd=lwd,...)
    mtext(ylab[1],side=3,at=1,cex=cex.lab)
    axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend),cex.axis=cex.axis)
    axis(2,at=c(0,0.5,1),labels=c("","",""),cex.axis=cex.axis)   
    
    abline(v=xlp)    
    for(i in 1:nx)
    {
      if(numeq(max(optmatrix[,i]),min(optmatrix[,i]))==TRUE)
      {
        mtext(c(round(min(optmatrix[,i])*0.75,digits=2+countdig(min(optmatrix[,i])*0.5)),round(1*(max(optmatrix[,i])),digits=2+countdig(1*(max(optmatrix[,i])))),round(1.25*max(optmatrix[,i]),digits=2+countdig(1.25*max(optmatrix[,i])))),side=2,at=c(0-(i-1)*(cex.axis*3/4)*1.25/26,0.5-(i-1)*(cex.axis*3/4)*1.25/26,1-(i-1)*(cex.axis*3/4)*1.25/26),col=col[i],cex=cex.axis,line=0.5)    
      }
      else
      {
        mtext(c(round(min(optmatrix[,i]),digits=2+countdig(min(optmatrix[,i]))),round(0.5*(max(optmatrix[,i])+min(optmatrix[,i])),digits=2+countdig(0.5*(max(optmatrix[,i])+min(optmatrix[,i])))),round(max(optmatrix[,i]),digits=2+countdig(max(optmatrix[,i])))),side=2,at=c(0-(i-1)*(cex.axis*3/4)*1.25/26,0.5-(i-1)*(cex.axis*3/4)*1.25/26,1-(i-1)*(cex.axis*3/4)*1.25/26),col=col[i],cex=cex.axis,line=0.5)
      }
      for(j in 1:length(xlp))
      {
        points(xlp[j],optplot[xlp[j],i],col="black",cex=cex)
      }
    }
    legend("topright",dimnames(out$Parameters)[[2]][1:nx],col=col[1:nx],lty=lty[1:nx],bty="n",cex=cex.lab,lwd=lwd)
 
  
    # right Plot  
      
    ##### First step to get every response on the same scale:
    targetvaluespos<-NULL
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    reoptplusdev<-matrix(NaN,ncol=ny,nrow=numbW)
    reoptminusdev<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
      
      #deviation
    if(standard==TRUE)
    {
      for(i in 1:ny)
      {
        reoptplot[,i]<-reoptmatrix[,i]-min(reoptmatrix[,i]-deviation[,i],tau[i])
        reoptplusdev[,i]<-reoptmatrix[,i]+deviation[,i]-min(reoptmatrix[,i]-deviation[,i],tau[i])
        reoptminusdev[,i]<-reoptmatrix[,i]-deviation[,i]-min(reoptmatrix[,i]-deviation[,i],tau[i])
        targetvaluespos[i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.5,ifelse(tau[i]>max(reoptmatrix[,i]+deviation[,i]),1,(tau[i]-min(c(reoptmatrix[,i]-deviation[,i],tau[i])))/(max(reoptmatrix[,i]+deviation[,i])-min(c(reoptmatrix[,i]-deviation[,i],tau[i])))))
        for(j in 1:numbW)
        { 
          reoptplot[j,i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.5,ifelse(tau[i]>max(reoptmatrix[,i]+deviation[,i]),reoptplot[j,i]/(tau[i]-min(reoptmatrix[,i]-deviation[,i])),reoptplot[j,i]/(max(reoptmatrix[,i]+deviation[,i])-min(reoptmatrix[,i]-deviation[,i],tau[i]))))
          reoptplusdev[j,i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.65,ifelse(tau[i]>max(reoptmatrix[,i]+deviation[,i]),reoptplusdev[j,i]/(tau[i]-min(reoptmatrix[,i]-deviation[,i])),reoptplusdev[j,i]/(max(reoptmatrix[,i]+deviation[,i])-min(reoptmatrix[,i]-deviation[,i],tau[i]))))
          reoptminusdev[j,i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.35,ifelse(tau[i]>max(reoptmatrix[,i]+deviation[,i]),reoptminusdev[j,i]/(tau[i]-min(reoptmatrix[,i]-deviation[,i])),reoptminusdev[j,i]/(max(reoptmatrix[,i]+deviation[,i])-min(reoptmatrix[,i]-deviation[,i],tau[i]))))
        }   
      }
      par(mar=c(4,3,3,2),mgp=c(3,1,0),bty=bty,las=las)
      matplot(xaxis1,reoptplot,type="l",col=col[(nx+1):(nx+ny)],lty=lty[(nx+1):(nx+ny)],xlab=xlab[2],ylab="",ylim=c(0,1.25),xaxt="n",yaxt="n",cex.lab=cex.lab,lwd=lwd,...)
      mtext(ylab[2],side=3,at=1,cex=cex.lab)
      axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend),cex.axis=cex.axis)
      axis(2,at=c(0,0.5,1),labels=c("","",""),cex.axis=cex.axis)   
      for(j in 1:ny)
      {
        polygon(c(1:numbW,numbW:1,1),c(reoptminusdev[,j],reoptplusdev[length(reoptplusdev[,j]):1,j],reoptminusdev[1,j]),col=rgb(col2rgb(col[nx+j])[1]/255,col2rgb(col[nx+j])[2]/255,col2rgb(col[nx+j])[3]/255,alpha=0.1+j*0.2/ny),border=FALSE)
      }
      for(i in 1:ny)                                                                                                                 
      {
        mtext(c(round(min(c(reoptmatrix[,i]-deviation[,i],tau[i])),digits=2+countdig(min(c(reoptmatrix[,i]-deviation[,i],tau[i])))),round((max(c(reoptmatrix[,i]+deviation[,i],tau[i]))+min(c(reoptmatrix[,i]-deviation[,i],tau[i])))*0.5,digits=2+countdig((max(c(reoptmatrix[,i]+deviation[,i],tau[i]))+min(c(reoptmatrix[,i]-deviation[,i],tau[i])))*0.5)),round(max(c(reoptmatrix[,i]+deviation[,i],tau[i])),digits=2+countdig(max(c(reoptmatrix[,i]+deviation[,i],tau[i]))))),side=2,at=c(0-(i-1)*(cex.axis*3/4)*1.25/26,0.5-(i-1)*(cex.axis*3/4)*1.25/26,1-(i-1)*(cex.axis*3/4)*1.25/26),col=col[nx+i],lty=lty[nx+i],cex=cex.axis,line=0.5)  
      } 
    }
    if(standard==FALSE)
    {
      for(i in 1:ny)
      {
        reoptplot[,i]<-reoptmatrix[,i]-min(reoptmatrix[,i],tau[i])
        targetvaluespos[i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.5,ifelse(tau[i]>max(reoptmatrix[,i]),1,(tau[i]-min(c(reoptmatrix[,i],tau[i])))/(max(reoptmatrix[,i])-min(c(reoptmatrix[,i],tau[i])))))
        for(j in 1:numbW)
        { 
          reoptplot[j,i]<-ifelse(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE,0.5,ifelse(tau[i]>max(reoptmatrix[,i]),reoptplot[j,i]/(tau[i]-min(reoptmatrix[,i])),reoptplot[j,i]/(max(reoptmatrix[,i])-min(reoptmatrix[,i],tau[i]))))
        }  
      }
      par(mar=c(4,3,3,2),mgp=c(3,1,0),bty=bty,las=las)
      matplot(xaxis1,reoptplot,type="l",col=col[(nx+1):(nx+ny)],lty=lty[(nx+1):(nx+ny)],xlab=xlab[2],ylab="",ylim=c(0,1.25),xaxt="n",yaxt="n",cex.lab=cex.lab,lwd=lwd,...)
      mtext(ylab[2],side=3,at=1,cex=cex.lab)
      axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend),cex.axis=cex.axis)
      axis(2,at=c(0,0.5,1),labels=c("","",""),cex.axis=cex.axis)   
        
      for(i in 1:ny)
      {
        if(numeq(tau[i],min(reoptmatrix[,i]))==TRUE&&numeq(tau[i],max(reoptmatrix[,i]))==TRUE)
        {
          mtext(c(round(ifelse(numeq(tau[i],0)==TRUE,-0.5,tau[i]-0.5*tau[i])),ifelse(numeq(tau[i],0)==TRUE,0,tau[i]),ifelse(numeq(tau[i],0)==TRUE,0.5,tau[i]+0.5*tau[i])),side=2,at=c(0-(i-1)*(cex.axis*3/4)*1.25/26,0.5-(i-1)*(cex.axis*3/4)*1.25/26,1-(i-1)*(cex.axis*3/4)*1.25/26),col=col[nx+i],cex=cex.axis,line=0.5)  
        }
        else
        {
          mtext(c(round(min(reoptmatrix[,i],tau[i]),digits=2+countdig(min(reoptmatrix[,i],tau[i]))),round(0.5*(max(reoptmatrix[,i],tau[i])+min(reoptmatrix[,i],tau[i])),digits=2+countdig(0.5*(max(reoptmatrix[,i],tau[i])+min(reoptmatrix[,i],tau[i])))),round(max(reoptmatrix[,i],tau[i]),digits=2+countdig(max(reoptmatrix[,i],tau[i])))),side=2,at=c(0-(i-1)*(cex.axis*3/4)*1.25/26,0.5-(i-1)*(cex.axis*3/4)*1.25/26,1-(i-1)*(cex.axis*3/4)*1.25/26),col=col[nx+i],cex=cex.axis,line=0.5)
        }
      }    
    }
        ######################
        ### Target Values ####
        ######################
        
    zaehler<-vector("list",length(tau))
    for(i in 1:length(tau))
    {
      counter<-NULL
      index<-1
      for(j in 1:length(tau))
      { 
        if(numeq(targetvaluespos[i],targetvaluespos[j])==TRUE && i!=j)
        {
          counter[index]<-j
          index<-index+1
        }
      }
      zaehler[[i]]<-sort(c(i,counter))
    }
    for(i in 1:length(tau))
    {
      point1<-NULL
      point2<-NULL
      for(j in 1:length(zaehler[[i]]))
      {
        if(i==zaehler[[i]][j])
        {
          point1<-c(1+(j-1)*(numbW-1)/length(zaehler[[i]]),1+j*(numbW-1)/length(zaehler[[i]]))
          point2<-c(targetvaluespos[i],targetvaluespos[i])
          lines(point1,point2,col=col[nx+i],lty=lty[nx+i],lwd=lwd)
          mtext(round(tau[i],digits=2+countdig(tau[i])),side=4,at=targetvaluespos[i]-(j-1)*(cex.axis*3/4)*1.25/26,col=col[nx+i],cex=cex.axis,line=-0.2)
        }
      }
    }  
    
        ##################################
        ####### End: Target Values #######
        ##################################

    
    yloc1<-c(0,1.25)
    abline(v=xlp,lty="solid")
    nam<-dimnames(out$Responses)[[2]]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau1[i],sep="=")
    }
    for(i in 1:ny)
    {
      for(j in 1:length(xlp))
      {
        points(xlp[j],reoptplot[xlp[j],i],col="black",cex=cex)#cols[nx+i],cex=3)
      }
    }
    legend("topright",nam,col=col[(nx+1):(nx+ny)],lty=lty[(nx+1):(nx+ny)],bty="n",cex=cex.lab,lwd=lwd)
  
    opt<-list(optp,reoptp)
  
    names(opt)<-list("ChosenParameters","ChosenResponses")

    return(opt)
}


