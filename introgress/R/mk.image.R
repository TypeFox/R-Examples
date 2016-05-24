mk.image <-
function(introgress.data=NULL, loci.data=NULL, marker.order=NULL, hi.index=NULL,
                   ind.touse=NULL, loci.touse=NULL, ylab.image="Individuals",
                   main.image="", xlab.h="population 2 ancestry", col.image=NULL,
                   pdf=TRUE, out.file="image.pdf"){
  if (is.data.frame(hi.index)==TRUE) hi.index<-hi.index[,2]	
  if (!is.numeric(hi.index)){
    stop("Please provide a numeric vector of hybrid indexes (hi.index).")
  }
  if(is.null(introgress.data)==TRUE)
    stop("error, introgress.data was not supplied")
  if(is.null(loci.data)==TRUE)
    stop("error, loci.data was not supplied")
  ## break up introgress.data if it is supplied as a list (i.e. if it was generated from
  ## prepare.data
  if (is.list(introgress.data)==TRUE)
    count.matrix<-introgress.data[[2]]
  else
    count.matrix<-introgress.data
  if (pdf==TRUE){
    pdf(file=out.file, height=9.5,width=7)
  }
  
  ## test dimension of loci.data to see if map information was supplied
  if(dim(loci.data)[2]>2){
    ## count number of markers per lg
    lg<-sort(unique(as.numeric(loci.data[,3])))
    marker.n<-numeric(length(lg))
    for(i in 1:length(lg)){
      marker.n[i] <- sum(loci.data[,3]==lg[i])
    }
  }
  if(is.null(loci.touse)){
    loci.touse<-1:(dim(count.matrix)[1])
  }
  if(is.null(ind.touse)){
    ind.touse<-1:(dim(count.matrix)[2])
  }
  if(is.null(col.image)){
    col.image<-c("#A1D99B", "#41AB5D", "#005A32")
  }
  if(is.null(marker.order)){
    marker.order<-loci.touse
  }
  else if(is.character(marker.order) || is.factor(marker.order)){
    tmp<-numeric(length(marker.order))
    for(i in 1:length(marker.order)){
      tmp[i] <- which(loci.data == as.character(marker.order)[i])
    }
    marker.order<-tmp
  }
  

  nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), c(6,1), c(1,1), respect=FALSE)  
  par(mar=c(5,4,1,1))

  ## can subset count.matrix and h individuals
  image(x=loci.touse, y=ind.touse,
        z=count.matrix[marker.order, order(hi.index[ind.touse])], col=col.image,
        breaks=c(-0.1,0.9,1.9,2.1), xlab="", ylab=ylab.image,
        main=main.image, axes=FALSE)
  mtext("Markers", 1, line=4)
  axis(2)
  ##axis(4)
  if(dim(loci.data)[2]==2){
    axis(1, at=seq(0.5, length(loci.touse)+0.5, 1), labels=FALSE)
    axis(1, at=loci.touse, tick=FALSE,
         labels=loci.data[,1][marker.order],
         cex.axis=0.3, las=2, line=1)
  }
  else{
    ## we know that there are at least linkage group data 
    axis(1, at=loci.touse, tick=FALSE,
         labels=loci.data[,1][marker.order],
         cex.axis=0.3, las=2, line=1)
    axis(1, at=c(0, cumsum(marker.n)+0.5), labels=FALSE)
    abline(v=c(0, cumsum(marker.n)+0.5), col="lightgray", lwd=0.6)
    axis(1, at=diff(c(0, cumsum(marker.n)))/2 +
         c(0, cumsum(marker.n)[-length(marker.n)]+0.5),
         tick=FALSE, labels=lg, cex.axis=0.5, line=0)
  }

  box()

  ## -- second plot -----------------------------------
  par(mar=c(5,0,1,1))

  n.inds<-length(ind.touse)
  ## subset hi.hindex appropriately if individuals are dropped

  plot(sort(hi.index[ind.touse]), 1:n.inds, axes=FALSE, xlab="", ylab="",
       ylim=c(1 + 0.0355 * n.inds, n.inds - 0.0355 * n.inds),
       ## kludge, ratio may depend on point size
       cex=0.6, type="n", xlim=0:1)
  
  mtext("Proportion", 1, line=3, cex=0.6)
  mtext(xlab.h, 1, line=4, cex=0.6)
  abline(v=0.5, col="lightgray")
  lines(sort(hi.index[ind.touse]), 1:n.inds)
  axis(1, at=c(0,0.5,1), cex.axis=0.6)
  ## axis(2)
  box()

  if(pdf==TRUE){
    dev.off()
  }
}

