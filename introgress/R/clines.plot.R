clines.plot <-
function(cline.data=NULL, marker.order=NULL, rplots=3, cplots=3,
                      pdf=TRUE, out.file="clines.pdf",colors=c("#005A32","#41AB5D"),
                      quantiles=FALSE,lb.cd=rep(0.025,3),ub.cd=rep(0.975,3),
                      lb.dh=rep(0.025,2),ub.dh=rep(0.975,2),cd=c("AA","Aa","aa"),dh=c("A","a")){
  if (is.null(cline.data)==TRUE)
    stop("error, input file was not provided")
  if (is.list(cline.data)==FALSE)
    stop("error, cline.data input file not in list format")
  if (length(colors)!=2) stop ("colors must be a vector of length two.")
  if (quantiles==TRUE & is.null(cline.data$Quantiles)==TRUE) stop ("error, quantile data not found.")
  ##set up output file
  if (pdf==TRUE) pdf(file=paste(out.file))
  par(mfrow=c(rplots,cplots))
  par(pty="s")
  ##break up data object cline.data
  ind.count.matrix<-cline.data[[8]]
  hi.index<-cline.data[[9]]
  loci.data<-cline.data[[10]]
  n.loci<-dim(loci.data)[1]
  ## if map data is provided clines are plotted in map order,
  ## otherwise they are plotted sequentially based on the cline.data
  ## object, marker order can also be specified by including a numeric
  ## or character vector giving the marker order
  if(is.null(marker.order) ){
      marker.order<-1:dim(loci.data)[1]
  }
  else if(is.character(marker.order) || is.factor(marker.order)){
    tmp<-numeric(n.loci)
    for(i in 1:n.loci){
      tmp[i] <- which(loci.data == as.character(marker.order)[i])
    }
    marker.order<-tmp
  }
  ## get quantile data for plotting if quantiles=TRUE
  if(quantiles==TRUE){
    gen.quantiles<-character(n.loci)
    for (i in marker.order){
      if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
        temp<-character(3)
        for (j in 1:length(temp)){
          if (cline.data$Quantiles[i,j] > ub.cd[j]) temp[j]<-"+"
          else if (cline.data$Quantiles[i,j] < lb.cd[j]) temp[j]<-"-"
          else temp[j]<-" "
        }
      }
      else {
        temp<-character(2)
        for (j in 1:length(temp)){
          if (cline.data$Quantiles[i,j] > ub.dh[j]) temp[j]<-"+"
          else if (cline.data$Quantiles[i,j] < lb.dh[j]) temp[j]<-"-"
          else temp[j]<-" "
        }
      }
      if (length(temp)==3){
        gen.quantiles[i]<-paste(cd[1],temp[1],"  ",cd[2],temp[2],"  ",cd[3],temp[3],sep="")
      }
      else gen.quantiles[i]<-paste(dh[1],temp[1],"  ",dh[2],temp[2],sep="")
    }
  }
  
  ## loop through loci in marker.order making plots
  for(i in marker.order){
    plot(0:1,0:1,type="n",xlab="Hybrid index",ylab="Pr(genotype)",axes=FALSE)
    ## add + and - for quantiles if quantiles==TRUE
    if(quantiles==TRUE){
      title(main=gen.quantiles[i],adj=0,cex.main=0.85)
    }
    ## determine if locus is dominant, haploid, or codominant
    if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
      if (is.null(cline.data[[5]])==FALSE){
        polygon.data.matrix<-rbind(hi.index,cline.data[[5]][[2]][i,],
                                   cline.data[[5]][[1]][i,],cline.data[[6]][[2]][i,],
                                   cline.data[[6]][[1]][i,],cline.data[[2]][i,],
                                   cline.data[[3]][i,])[,order(hi.index)]
        polygon.data.matrix.rev<-rbind(hi.index,cline.data[[5]][[2]][i,],
                                       cline.data[[5]][[1]][i,],cline.data[[6]][[2]][i,],
                                       cline.data[[6]][[1]][i,],cline.data[[2]][i,],
                                       cline.data[[3]][i,])[,order(-hi.index)]
        polygon(c(polygon.data.matrix[1,],polygon.data.matrix.rev[1,]),
                c(polygon.data.matrix[2,],polygon.data.matrix.rev[3,]),
                col=colors[1], border=NA)
        polygon(c(polygon.data.matrix[1,],polygon.data.matrix.rev[1,]),
                c(polygon.data.matrix[4,],polygon.data.matrix.rev[5,]),
                col=colors[2],border=NA)
        lines(polygon.data.matrix[1,],polygon.data.matrix[6,],lty=1)
        lines(polygon.data.matrix[1,],polygon.data.matrix[7,],lty=2)
      }
      else {
        polygon.data.matrix<-rbind(hi.index,cline.data[[2]][i,],
                                   cline.data[[3]][i,])[,order(hi.index)]
        polygon.data.matrix.rev<-rbind(hi.index,cline.data[[2]][i,],
                                       cline.data[[3]][i,])[,order(-hi.index)]
        lines(polygon.data.matrix[1,],polygon.data.matrix[2,],lty=1)
        lines(polygon.data.matrix[1,],polygon.data.matrix[3,],lty=2)
      }
      points(hi.index,ind.count.matrix[i,]/2)
      text(0.75,0.95, loci.data[i,1], pos=1, cex=0.8)
      if (is.null(cline.data[[5]])==FALSE)
        text(0.75,0.85, paste("P=",signif(as.numeric(cline.data[[1]][i,4]),3)),
             pos=1, cex=0.8)
      axis(4, at=c(0,.5, 1), labels=hist(ind.count.matrix[i,], breaks=-1:2, plot=FALSE)$counts,
           las=1, cex.axis=0.6)
      axis(1, at=c(0,.5, 1))
      axis(2, at=c(0,.5, 1))
      box()
    }
    ## for dominant loci		
    else if (loci.data[i,2]=="D" | loci.data[i,2]=="d" | loci.data[i,2]=="H" | loci.data[i,2]=="h"){
      if (is.null(cline.data[[5]])==FALSE){
        polygon.data.matrix<-rbind(hi.index, cline.data[[5]][[2]][i,],
                                   cline.data[[5]][[1]][i,],cline.data[[2]][i,])[,order(hi.index)]
        polygon.data.matrix.rev<-rbind(hi.index, cline.data[[5]][[2]][i,],
                                       cline.data[[5]][[1]][i,],cline.data[[2]][i,])[,order(-hi.index)]
        polygon(c(polygon.data.matrix[1,],polygon.data.matrix.rev[1,]),
                c(polygon.data.matrix[2,],polygon.data.matrix.rev[3,]),
                col=colors[2], border=NA)
        lines(polygon.data.matrix[1,],polygon.data.matrix[4,],lty=1)
      }
      else {
        polygon.data.matrix<-rbind(hi.index,cline.data[[2]][i,])[,order(hi.index)]
        polygon.data.matrix.rev<-rbind(hi.index,cline.data[[2]][i,])[,order(-hi.index)]
        lines(polygon.data.matrix[1,],polygon.data.matrix[2,],lty=1)
      }
      points(hi.index,ind.count.matrix[i,])
      text(0.75,0.95, loci.data[i,1], pos=1, cex=0.8)
      if (is.null(cline.data[[5]])==FALSE)
        text(0.75,0.85, paste("P=",signif(as.numeric(cline.data[[1]][i,4]),3)),
             pos=1, cex=0.8)
      axis(4, at=c(0, 1), labels=hist(ind.count.matrix[i,], breaks=c(-1,0.5,2), plot=FALSE)$counts,
           las=1, cex.axis=0.6)
      axis(1, at=c(0,.5, 1))
      axis(2, at=c(0,.5, 1))
      box()
    }		
  }
  if (pdf==TRUE) dev.off()
}

