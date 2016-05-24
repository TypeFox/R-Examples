sample.plot <- function(data, labels=NULL, plot=TRUE, ...){

    if (is.vector(data)){
        data2 <- matrix(data, nrow=1, byrow=T)
        rownames(data2) <- ""
        colnames(data2) <- names(data)
        data <- data2
        rm(data2)
    }
    
  ##Have to add the name for apply function
  data2<-data.frame(matrix(NA, ncol=ncol(data)+1, nrow=nrow(data)))
  data2[,1]<-rownames(data)
  data2[,2:ncol(data2)]<-data
  colnames(data2)<-c("probe", colnames(data))
  data<-data2
  rm(data2)

  ##For each variable
  out<-apply(data,1,exp_grp<-function(x){
    
    id<-x[1]
    x<-as.numeric(x[-1])
    s.name<-colnames(data)[-1]
    if (!is.null(labels)){
    
      x.g<-x.name<-x.col<-x.leg<-c()
      col<-rainbow(n=length(unique(labels)))
      cpt<-1

      for (i in unique(labels)){
        tmp<-x[which(labels==i)]
        x.g<-c(x.g,tmp)
        x.name<-c(x.name, s.name[which(labels==i)])
        x.col<-c(x.col, rep(col[cpt],length(tmp)))
        x.leg<-c(x.leg,i)
        cpt<-cpt+1
      }

      if (plot){
        dev.new()
        barplot(x.g,names.arg=x.name,axis.lty=1,las=2,col=x.col, ylab='Expression Level',main = paste("Expression level of",id), cex.axis=0.7, cex.lab=0.7, cex.main=0.8, cex.names=0.6, ylim=c(0,(max(x,na.rm=TRUE)+1)),...)
        legend(x=-1,y=(max(x,na.rm=TRUE)+1),x.leg,fill=col, bty="n", cex=0.7, text.col="gray50")
      }
    }
    else{
      if (plot){
        dev.new()
        barplot(x,names.arg=s.name, axisnames=TRUE,las=2,col="blue", ylab='Expression Level',main = paste("Expression level of",id), cex.axis=0.7, cex.lab=0.7, cex.main=0.8, cex.names=0.6, ...)
      }
    }
  })

  return (out)
}

distrib.plot <- function(data, labels=NULL, plot=TRUE, ...){
    if (is.vector(data)){
        data <- matrix(data, nrow=1)
    }
    
    ##Have to add the name for apply function
    data2<-data.frame(matrix(NA, ncol=ncol(data)+1, nrow=nrow(data)))
    data2[,1]<-rownames(data)
    data2[,2:ncol(data2)]<-data
    colnames(data2)<-c("probe", colnames(data))
    data<-data2
    rm(data2)
    
    ##For each variable
    apply(data,1,exp_grp<-function(x){
        dev.new()
        id<-x[1]
        x<-as.numeric(x[-1])
        
        if (!is.null(labels)){
            col<-rainbow(n=length(unique(labels)))
            cpt<-1
            
            maxcounts<-NULL
            for (i in unique(labels)){
                x.tmp<-x[which(labels==i)]
                if(cpt==1){
                    plot(density(x.tmp, na.rm=TRUE), col=col[cpt], type="l", lwd=2, cex.axis=0.7, cex.lab=0.7, cex.main=0.8,main = paste("Expression level of",id),xlab='Expression Level',ylab='Density', ...)
                }
                else{
                    points(density(x.tmp, na.rm=TRUE), col=col[cpt], type="l", lwd=2)
                }
                
                cpt<-cpt+1
            }
            legend(x=min(x),y=maxcounts,unique(labels),fill=col, bty="n", cex=0.7, text.col="gray50", ...)
        }
        else{
            plot(density(x, na.rm=TRUE), cex.lab=0.7, col="blue",cex.axis=0.7,cex.main=0.8,main = paste("Expression level of",id),xlab='Expression Level',ylab='Density', ...)
        }
    })
}

probePlots <-  function(abatch, path, pbsList, labAxisProbes=TRUE, labAxisArrays=TRUE, legendArrays=TRUE, legendProbes=TRUE, cex.axis=0.9, cex.legend=0.8, pdfName){

    ## check if arguments are ok
    if(missing(abatch) && missing(path))
        stop("** ERROR : 'abatch' or 'path' arguments are empty")
    
    if(missing(pbsList))
        stop("** ERROR : 'pbsList' argument is missing")
    if(class(pbsList)!="character")
        stop("** ERROR : 'pbsList' argument is not a character class")
    
    if(!missing(abatch)){
        if(class(abatch) ==  "AffyBatch")
            abatch <- abatch
        else
            stop("** ERROR : 'abatch' argument is not an AffyBatch class")
    }
    
    if(!missing(path)){
        if(class(path) ==  "character")
            abatch <- ReadAffy(celfile.path=path)
        else
            stop("** ERROR : 'path' argument is not a character class")
    }
    
    if(!missing(pdfName)){
        if(class(pdfName) ==  "character")
            pdf(pdfName, width=9, height=8)
        else
            stop("** ERROR : 'pdfName' argument is not a character")
    }
    
    
    ## select the probesets in the affybatch
    pbs <- probeset(abatch, pbsList)
    
    ## prepare the plot
    par(mfrow=c(2,3), bg = "gray80", mai=c(5, 4, 4, 2), mar=c(5, 4, 4, 2), oma=c(0,0,2,0))
    
    for(i in 1:length(pbs)){
        pmi <- pm(pbs[[i]])
        mmi <- mm(pbs[[i]])
        
        ## INTER CHIPS PLOT
        plot(pbs[[i]], type="l", main=paste("Probeset ", names(pbs)[i], "\nInter arrays\nxaxis = probes", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(pmi)), lty=1:5)
        if(!is.null(rownames(pmi)) && labAxisProbes){
            xlabels <- rownames(pmi)
            axis(1, 1:nrow(pmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(is.null(rownames(pmi)) && labAxisProbes){
            xlabels <- as.character(1:nrow(pmi))
            axis(1, 1:nrow(pmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(!labAxisProbes)
            axis(1, 1:nrow(pmi), labels=FALSE, las=2, cex.axis=cex.axis)
        axis(2)
        box()
        if(!is.null(colnames(pmi)) && legendArrays)
            legend("topright", colnames(pmi), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
        
        ## INTER PROBES PLOT - PM
        matplot(t(pmi), type="l", main=paste("Probeset ", names(pbs)[i], "\nInter probes - Perfect Match\nxaxis = arrays", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(pmi)), lty=1:5)
        if(!is.null(colnames(pmi)) && labAxisArrays){
            xlabels <- colnames(pmi)
            axis(1, 1:ncol(pmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(is.null(colnames(pmi)) && labAxisArrays){
            xlabels <- as.character(1:ncol(pmi))
            axis(1, 1:ncol(pmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(!labAxisArrays)
            axis(1, 1:ncol(pmi), labels=FALSE, las=2, cex.axis=cex.axis)
        axis(2)
        box()
        if(!is.null(rownames(pmi)) && legendProbes)
            legend("topright", rownames(pmi), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
        if(is.null(rownames(pmi)) && legendProbes)
            legend("topright", as.character(1:nrow(pmi)), lty=1:5, col=rainbow(nrow(pmi)), cex=cex.legend)
        
        ## INTER PROBES PLOT - MM
        matplot(t(mmi), type="l", main=paste("Probeset ", names(pbs)[i], "\nInter probes - Mis Match\nxaxis = arrays", sep=""), axes=FALSE, ylab="Intensity", xlab="", col=rainbow(nrow(mmi)), lty=1:5)
        if(!is.null(colnames(mmi)) && labAxisArrays){
            xlabels <- colnames(mmi)
            axis(1, 1:ncol(mmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(is.null(colnames(mmi)) && labAxisArrays){
            xlabels <- as.character(1:ncol(mmi))
            axis(1, 1:ncol(mmi), xlabels, las=2, cex.axis=cex.axis)
        }
        if(!labAxisArrays)
            axis(1, 1:ncol(mmi), labels=FALSE, las=2, cex.axis=cex.axis)
        axis(2)
        box()
        if(!is.null(rownames(mmi)) && legendProbes)
            legend("topright", rownames(mmi), lty=1:5, col=rainbow(nrow(mmi)), cex=cex.legend)
        if(is.null(rownames(mmi)) && legendProbes)
            legend("topright", as.character(1:nrow(mmi)), lty=1:5, col=rainbow(nrow(mmi)), cex=cex.legend)
        
    }
    if(!missing(pdfName)){
        if(class(pdfName) ==  "character")
            dev.off()
    }
}


genes.selection <- function(data, thres.diff, thres.num, probs=0.25){
    ## check the parameters
    if (missing(thres.diff) && missing(thres.num)){
        stop("** Stop. No method found.")
    }
    if (!missing(thres.diff) && !missing(thres.num)){
        stop("** Stop. Choose one of the two options - thres.diff or thres.num")
    }
    
    diff.quantile <- function(x, probs=1/length(x)){
        vv <- quantile(x, probs=c(probs, 1-probs), na.rm = TRUE)
        return(vv[2] - vv[1])
    }
    
    rangeValues<-apply(data, 1, diff.quantile, probs=probs)
    
    if (!missing(thres.diff)){
        ind<-which(rangeValues>=thres.diff)
        genesList<-rownames(data[ind,])    
    }
    else if (!missing(thres.num)){
        genesList<-names(sort(rangeValues, decreasing=TRUE)[1:thres.num])
    }
    return (genesList)
}

