#
# power.plot.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified October, 2013
# first written January, 2013
# Contains: power.plot.R
#

# cross.saturate
#
# DESCRIPTION:
#  Saturate an existing genetic map by adding markers derived from expression
# OUTPUT:
#  An object of class cross
#
power.plot <- function(cross1,cross2,scores,qtlThr=5,nPheno=500,verbose=FALSE,...){
  if(missing(scores)){
    if(missing(cross1)){
      stop("No object of class cross provided (cross1).")
    }
    if(missing(cross2)){
      stop("No object of class cross provided (cross2).")
    }
    if (!any(class(cross1) == "cross")){
      stop("Input should have class \"cross\" (cross1).")
    }
    if (!any(class(cross2) == "cross")){
      stop("Input should have class \"cross\" (cross2).")
    }
    if(!is.numeric(qtlThr) || qtlThr<0){
      stop("qtlThr should be a numeric value bigger than 0.")
    }
    if(nphe(cross1)!=nphe(cross2)){
       stop("Difference in phenotype number between crosses.")
    }
    if(!is.numeric(nPheno) || nPheno<2 || nPheno>nphe(cross1)){
      stop("nPheno should be a numeric value between 2 and nphe(cross1).")
    }
    if(is.null(cross1$geno[[1]]$prob)){
      cross1 <- calc.genoprob(cross1)
    }
    if(is.null(cross2$geno[[1]]$prob)){
      cross2 <- calc.genoprob(cross2)
    }
    markers        <- sort(sample(1:nphe(cross1),nPheno))
    qtlScores1     <- scanone(cross1,pheno.col=markers,verbose=verbose,...)
    qtlScores2     <- scanone(cross2,pheno.col=markers,verbose=verbose,...)
        
  }else if(class(scores)=="scores"){
    qtlScores1     <- scores[[1]]
    qtlScores2     <- scores[[2]]
  }else{
    stop("This is not a valid scores object!")
  }
  maxes1   <- apply(qtlScores1[,3:ncol(qtlScores1)],2,max)
  maxes2   <- apply(qtlScores2[,3:ncol(qtlScores2)],2,max)
  if(verbose){
    vals           <- cbind(markers,maxes1,maxes2)
    colnames(vals) <- c("#phenotype","max LOD in cross1","max LOD in cross2")
    print(vals)
  }
  colorNr      <- as.numeric(maxes1>maxes2)+1
  colorCols    <- c("blue","red")
  
  plot(maxes1,maxes2,pch=20,cex=0.3,xlab="LOD scores on the original map",ylab="LOD scores on the saturated map",col=colorCols[colorNr],log="xy")
  abline(qtlThr-qtlThr,1,col="red") #diagonal line
  abline(v=qtlThr,lty=2,col="grey")   #dotted line showing threshold used for x axis
  abline(h=qtlThr,lty=2,col="grey")  #dotted line showing threshold used for y axis
  
  significant  <- which(maxes1>qtlThr)  #nr of significant QTLs in cross 1
  significant2 <- which(maxes2>qtlThr) #nr of significant QTLs in cross 2
  stronger     <- sum(maxes1[significant]<(maxes2[significant])) #nr of QTLs showubg increase in power
  
  #verbose output
  cat("\n##################################\n")
  cat("Number significant QTLs:\n\t- original cross:",length(significant),"\n\t- new cross:     ",length(significant2),"\n\n")
  cat(stronger,"which is:",stronger/length(significant),"% of significant qtls show increase in power\n")
  cat(length(significant2)-length(significant),"new QTLs found\n") # nr of gained QTLs
  cat("\n##################################\n")
  
  scores <- list(qtlScores1,qtlScores2)
  class(scores) <- "scores"
  invisible(scores)
}

