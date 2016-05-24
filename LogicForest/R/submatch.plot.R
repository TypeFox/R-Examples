submatch.plot <-
function(fit, pred.nms, pis, preds, size, color=FALSE)
{  
  allPIs<-unique(names(fit$PI.frequency))
  if(missing(pred.nms)) {pred.nms<-paste("X", 1:preds, sep="")}
  matches<-persist.match(fit=fit, pred.nms=pred.nms, preds=preds)
  unq.sz<-matches$uniq.sz#unique sizes that were matched
  all.sz<-matches$all.sz#sizes observed in forest
  num.sz<-matches$numsizes 
  mat<-matches$num.matches
  colr<-colors()
  if (color==FALSE) {colvec<-c(271,321,281,331,291,341,301,351,311,360)}
  else {colvec<-c(552,132,254,550,149,74,547,100,121,393,498,142,455,31,102,116,115,638,598,657)}
  if (length(size)>1)
    {
    lo<-layout(matrix(c(2,0,5,7,1,3,4,6), 2, 4, byrow=TRUE), c(3,1,3,1), c(1,3), TRUE) 
    for (i in 1:length(size))
      { 
      sz<-which(unq.sz%in%size[i])
      mch.mat<-mat[[sz]]
      loc<-min(which(all.sz>=size[i]))
      if (length(loc)==0) {ellig.PIs<-sum(num.sz)}
      else {ellig.PIs<-sum(num.sz[loc:length(num.sz)])}
      cols<-ncol(mch.mat)
      pi<-min(pis, nrow(mch.mat))
      top.mat<-order(rowSums(mch.mat), decreasing=TRUE)[1:pi]
      plotmat<-mch.mat[c(top.mat),]
      rnms<-rownames(plotmat) 
      primfreq.match<-fit$PI.frequency[which(names(fit$PI.frequency)%in%rnms)]
      if (is.vector(plotmat)) {PIfreq<-sum(primfreq.match); SZfreq<-as.vector(primfreq.match)}
      else {
        for(k in 1:length(primfreq.match))
          {
          id<-which(rownames(plotmat)==names(primfreq.match[k]))
          plotmat[id,1]<-primfreq.match[k]
          }
        PIfreq<-rowSums(plotmat);  SZfreq<-colSums(plotmat)}
        top1<-max(PIfreq)
        top2<-max(SZfreq)
        tit<-paste("Persistance for size", size[i], sep=" ")
        par(mar=c(3,3,1,1), cex.axis=0.95, las=1, mai=c(1, 1, 0.25,0.25))
        if (is.vector(plotmat)) 
        {
        x<-c(1:cols)
        y<-rep(1, cols)
        symbols(x, y, circles=plotmat/ellig.PIs, ylab="",xlab="", main=tit, inches=0.2, xlim=c(0,(cols+.5)), 
        ylim=c(0,(pi+.5)), xaxt="n", yaxt="n", bg=colr[colvec[1]], fg=colr[colvec[1]])
        names(plotmat)<-names(num.sz)[loc:max(unq.sz)]
        abline(h=1, col=colr[colvec[1]], lty=3)
        axis(1, at=c(1:cols), labels=names(plotmat), cex.axis=0.75)
        axis(2, at=1, labels=rownames(mch.mat)[c(top.mat),], cex.axis=0.75)
        }
      else 
        {
        rownames(plotmat)<-rownames(mch.mat)[c(top.mat)]
        circs<-as.vector(t(plotmat))
        x<-rep(1:cols, pi)
        y<-rep(pi:1, each=cols, times=1)
        symbols(x, y, circles=circs/ellig.PIs, ylab="",xlab="", main=tit, inches=.2,xlim=c(0,(cols+.5)), 
               ylim=c(0,(pi+.5)), xaxt="n", yaxt="n", bg=rep(colr[colvec[1:pi]], each=cols, times=1),
               fg=rep(colr[colvec[1:pi]], each=cols, times=1))
        abline(h=c(pi:1), col=colr[colvec[1:pi]], lty=3)
        axis(1, at=c(1:cols), labels=colnames(plotmat))
        axis(2, at=c(pi:1), labels=rownames(plotmat))
        }  
      par(mar=c(0,3,1,1), mai=c(.1, 1, .1, .1))
      barplot(SZfreq, axes=FALSE, ylim=c(0,top2), space=0, axisnames=FALSE)
      par(mar=c(3,0,1,1), mai=c(1, .1, .1, .1))
      barplot(sort(PIfreq, decreasing=FALSE), axes=FALSE, xlim=c(0,top1), space=0, axisnames=FALSE, horiz=TRUE)
      }
    }
  else 
    {
    sz<-which(unq.sz%in%size)
    all.sz<-matches$all.sz
    num.sz<-matches$numsizes
    mat<-matches$num.matches[[sz]]
    rnms<-rownames(mat)
    primfreq.match<-fit$PI.frequency[which(names(fit$PI.frequency)%in%rnms)]
    for(k in 1:length(primfreq.match))
      {
      id<-which(rownames(mat)==names(primfreq.match[k]))
      mat[id,1]<-primfreq.match[k]
      }
    num.sz[size]<-colSums(mat)[1]
    loc<-min(which(all.sz>=size))
    if (length(loc)==0) {ellig.PIs<-sum(num.sz)}
    else {ellig.PIs<-sum(num.sz[loc:length(num.sz)])}
    cols<-ncol(mat)
    pi<-min(pis, nrow(mat))
    top.mat<-order(rowSums(mat), decreasing=TRUE)[1:pi]
    plotmat<-mat[c(top.mat),]
    PIfreq<-rowSums(plotmat)
    SZfreq<-colSums(plotmat)
    top1<-max(PIfreq)
    top2<-max(SZfreq)
    lo<-layout(matrix(c(2,0,1,3), 2, 2, byrow=TRUE), c(3,1), c(1,3), TRUE) 
    tit<-paste("Persistance for size", size, sep=" ")
    par(mar=c(3,3,1,1), cex.axis=0.75, las=1, mai=c(1, 1, 0.25,0.25))
    if (is.vector(plotmat)) 
      {
      x<-c(1:cols)
      y<-rep(1, cols)
      symbols(x, y, circles=plotmat/ellig.PIs, ylab="",xlab="", main=tit, inches=0.2, xlim=c(0,(cols+.5)), 
              ylim=c(0,(pi+.5)), xaxt="n", yaxt="n", bg=colr[colvec[1]], fg=colr[colvec[1]])
      names(plotmat)<-names(num.sz)[loc:length(num.sz)]
      abline(h=1, col=colr[colvec[1]], lty=3)
      axis(1, at=c(1:cols), labels=names(plotmat), cex.axis=.85)
      axis(2, at=1, labels=rownames(mat)[c(top.mat),], cex.axis=.85)
      }
    else 
      {
      rownames(plotmat)<-rownames(mat)[c(top.mat)]
      circs<-as.vector(t(plotmat))
      x<-rep(1:cols, pi)
      y<-rep(pi:1, each=cols, times=1)
      symbols(x, y, circles=circs/ellig.PIs, ylab="",xlab="", main=tit, inches=.2,xlim=c(0,(cols+.5)), 
      ylim=c(0,(pi+.5)), xaxt="n", yaxt="n", bg=rep(colr[colvec[1:pi]], each=cols, times=1),
      fg=rep(colr[colvec[1:pi]], each=cols, times=1))
      abline(h=c(pi:1), col=colr[colvec[1:pi]], lty=3)
      axis(1, at=c(1:cols), labels=colnames(plotmat), cex.axis=.85)
      axis(2, at=c(pi:1), labels=rownames(plotmat), cex.axis=.85)
      }
    top<-max(num.sz)
    par(mar=c(0,3,1,1), mai=c(.1, 1, .1, .1))
    if (min(all.sz)==1) {
      barplot(num.sz[sz:max(all.sz)], axes=FALSE, ylim=c(0,top), space=0, axisnames=FALSE, add=FALSE, col="grey90")
      par(mar=c(0,3,1,1), mai=c(.1, 1, .1, .1))
      barplot(SZfreq, axes=FALSE, ylim=c(0,top), space=0, axisnames=FALSE, add=TRUE, col="grey40")
      }
    if (min(all.sz)!=1){
      barplot(num.sz[(sz-1):(max(all.sz)-1)], axes=FALSE, ylim=c(0,top), space=0, axisnames=FALSE, add=FALSE, col="grey90")
      par(mar=c(0,3,1,1), mai=c(.1, 1, .1, .1))
      barplot(SZfreq, axes=FALSE, ylim=c(0,top), space=0, axisnames=FALSE, add=TRUE, col="grey40")
      }
    par(mar=c(3,0,1,1), mai=c(1, .1, .1, .1))
    PIfreqs<-rowSums(mat)
    col1<-rep("grey40", pis)
    if(length(PIfreqs)>pis) {col2<-rep("grey90", length(PIfreqs)-pis); colr<-c(col2, col1)}
    else {colr<-col1}
    topsz<-max(PIfreqs)
    barplot(sort(PIfreqs, decreasing=FALSE), axes=FALSE, xlim=c(0, topsz), col=colr, space=0, axisnames=FALSE, horiz=TRUE)
  }
}
