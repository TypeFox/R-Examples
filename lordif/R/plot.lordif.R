plot.lordif <-
function(x,labels=c("Reference","Focal"),width=7,height=7,...) {
    ndif<-sum(x$flag)
    if (ndif==0) stop(paste(deparse(substitute(x))), " contains no items flagged for DIF")
    if (ndif==x$ni) stop("all items in ", paste(deparse(substitute(x))), " have been flagged for DIF")
    if (x$ng != length(labels)) labels<-paste("Group",1:x$ng)
    sumpp<-function(pp) {
      ws<-rowSums(pp*(col(pp)-1))
      return(ws)
    }
    maxcat<-ncol(x$ipar.sparse)
    sysname<-Sys.info()[["sysname"]]
    if(sysname=="Windows") {
      dev.new(width=width,height=height,record=TRUE)
    } else if (sysname=="Linux") {
      dev.new(width=width,height=height)
      par(ask=TRUE)
    } else {
      dev.new(width=width,height=height)
    }
    par(mfrow=c(1,1))
    theta<-seq(x$options$minTheta,x$options$maxTheta,x$options$inc)
    difitems<-(1:x$ni)[x$flag]
    difselections<-x$selection[x$flag]
    itemnames<-row.names(x$ipar.sparse)
    gpar<-array(NA,c(ndif,maxcat,x$ng))
    cpar<-as.matrix(x$ipar.sparse[1:(x$ni-ndif),])
    pp<-array(NA,c(length(theta),ndif,maxcat,x$ng))
    gtheta<-split(x$calib.sparse$theta,x$group)
    gdensity<-matrix(0,length(theta),x$ng)
    for (i in 1:x$ng) {
      gdensity[,i]<-density(unlist(gtheta[names(table(x$group))[i]]),n=length(theta),from=x$options$minTheta,to=x$options$maxTheta,bw=.25)$y
    }
    plot(theta,gdensity[,1],type="l",xlab="theta",ylab="Density",ylim=c(0,max(gdensity)),lty=1,col=1,main="Trait Distributions",...)
    for (g in 2:x$ng) {
      lines(theta,gdensity[,g],lty=g,col=g)
    }
    legend("topright",labels,lty=1:x$ng,col=1:x$ng,bg="white")
    par(mfrow=c(2,2))
    for (i in 1:length(difitems)) {
      ncat<-x$ncat[difitems[i]]
      plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Item Score",main=paste0("Item True Score Functions - Item ",difselections[i]),...)
      for (g in 1:x$ng) {
        gpar[i,,g]<-unlist(x$ipar.sparse[which(itemnames==paste0("I",difselections[i],".",g)),])
        if(x$options$model=="GPCM") pp[,i,1:ncat,g]<-probgpcm(theta,gpar[i,1,g],gpar[i,2:ncat,g]) else pp[,i,1:ncat,g]<-probgrm(theta,gpar[i,1,g],gpar[i,2:ncat,g])
        lines(theta,sumpp(pp[,i,1:ncat,g]),lty=g,col=g)
      }
      legend("bottomright",labels,lty=1:x$ng,col=1:x$ng,cex=0.7,bg="white")
      chi12<-paste(x$stats[difitems[i],"df12"],")=",x$stats[difitems[i],"chi12"],sep="")
      pseudo12<-x$stats[difitems[i],paste("pseudo12.",x$options$pseudo.R2,sep="")]
      beta12<-round(x$stats[difitems[i],"beta12"],4)
      chi13<-paste(x$stats[difitems[i],"df13"],")=",x$stats[difitems[i],"chi13"],sep="")
      pseudo13<-x$stats[difitems[i],paste("pseudo13.",x$options$pseudo.R2,sep="")]
      chi23<-paste(x$stats[difitems[i],"df23"],")=",x$stats[difitems[i],"chi23"],sep="")
      pseudo23<-x$stats[difitems[i],paste("pseudo23.",x$options$pseudo.R2,sep="")]
      text(min(theta),ncat-1,substitute(paste("Pr(",chi[12]^2,",",chi12,",",R[12]^2,"=",pseudo12,",",Delta,"(",beta[1],")=",beta12,sep="")),adj=c(0,1),cex=0.8)
      text(min(theta),(ncat-1)*.9,substitute(paste("Pr(",chi[13]^2,",",chi13,",",R[13]^2,"=",pseudo13,sep="")),adj=c(0,1),cex=0.8)
      text(min(theta),(ncat-1)*.8,substitute(paste("Pr(",chi[23]^2,",",chi23,",",R[23]^2,"=",pseudo23,sep="")),adj=c(0,1),cex=0.8)
      plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Item Score",main="Differences in Item True Score Functions",...)
      for (g in 2:x$ng) {
        lines(theta,abs(sumpp(pp[,i,1:ncat,1])-sumpp(pp[,i,1:ncat,g])),lty=g,col=g)
      }
      plot(theta,seq(0,1,along.with=theta),type="n",xlab="theta",ylab="Probability",main="Item Response Functions",...)
      for (g in 1:x$ng) {
        for (k in 1:ncat) {
          lines(theta,pp[,i,k,g],lty=g,cex=0.1,col=g)
        }
      }
      for (g in 1:x$ng) {
        text(x$options$minTheta,.8-(g-1)*par()$cxy[2],paste(round(gpar[i,,g][!is.na(gpar[i,,g])],2),collapse=", "),col=g,adj=c(0,0),cex=0.8)
        for (k in 2:ncat) {
          if (!is.na(gpar[i,k,g])) text(gpar[i,k,g],0,"|",col=g)
        }
      }
      plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Size",main="Impact (Weighted by Density)",...)
      for (g in 2:x$ng) {
        lines(theta,gdensity[,g]*abs(sumpp(pp[,i,1:ncat,1])-sumpp(pp[,i,1:ncat,g])),lty=g,col=g)
      }
    }
    par(mfrow=c(1,2))
    plot(theta,seq(0,sum(!is.na(x$ipar))-x$ni,along=theta),xlab="theta",ylab="TCC",type="n",main="All Items",...)
    for (g in 1:x$ng) {
      apar<-rbind(cpar,gpar[,,g])
      lines(theta,tcc(apar[,1],apar[,-1,drop=F],theta,model=x$options$model),lty=g,col=g)
    }
    legend("bottomright",labels,lty=1:x$ng,col=1:x$ng,bg="white")
    plot(theta,seq(0,sum(!is.na(gpar[,,1]))-ndif,along=theta),xlab="theta",ylab="TCC",type="n",main="DIF Items",...)
    for (g in 1:x$ng) {
      lines(theta,tcc(gpar[,1,g],matrix(gpar[,-1,g],nrow=ndif),theta,model=x$options$model),lty=g,col=g)
    }
    legend("bottomright",labels,lty=1:x$ng,col=1:x$ng,bg="white")
    layout(matrix(c(1,2),ncol=2),widths=c(1,2))
    boxplot(x$calib$theta-x$calib.sparse$theta,col = "light grey")
    difference<-x$calib$theta-x$calib.sparse$theta
    plot(x$calib$theta,difference,type="n",xlab="initial theta",ylab="initial - purified",...)
    abline(h=0)
    abline(h=mean(x$calib$theta-x$calib.sparse$theta),lty=2)
    for (i in 1:x$ng) {
      points(x$calib$theta[x$group==as.numeric(names(table(x$group))[i])],difference[x$group==as.numeric(names(table(x$group))[i])],col=i,pch=i)
    }
    legend("topright",labels,pch=1:x$ng,col=1:x$ng,bg="white")
  }
