colorhcplot <-
function (hc, fac, hang = 0.1, main = "Cluster Dendrogram", color = T, lab.cex = 0.9, ax.low="auto", ax.top="auto", lwd=3, window.adj=1,tick.at=5) {
  
  #reorder leaves (samples) based on hc$order 
  #merges matrix (mgs) will include a third column -> middle point of each segment
  m_ord <- -1*hc$order
  mgs <- as.vector(hc$merge)
  for (i in 1:length(mgs)) {
    if (mgs[i] <0) { mgs[i] <- as.numeric(-1* which(m_ord == mgs[i]))}
  }
  mgs <- matrix(mgs, ncol=2)
  mgs <- cbind(mgs, rep(NA,nrow(mgs)))
  for (step in 1:nrow(mgs)) {
    if(mgs[step,1]<0 & mgs[step,2]<0) {
      mgs[step,3] <- (-1)*mean(mgs[step,1:2])
    } else if (mgs[step,1] * mgs[step,2]< 0) {
      min <- which(mgs[step,1:2]<0)
      plu <- which(mgs[step,1:2]>0)
      mgs[step,3] <- ((-1*mgs[step,min])+mgs[mgs[step,plu],3])/2
    } else {
      mgs[step,3] <- (mgs[mgs[step,1],3]+mgs[mgs[step,2],3])/2
    }
  }
  
  #cluster lookup and cbind in merges table (mgs)
  #mgs matrix will include two additional cols -> color/level of each
  #element involved in the merge. If no match, value is 8=gray.
  mgs<-mgs[,1:3] #just in case, remove unwanted columns
  mgs <- cbind(mgs, matrix(NA,nrow=nrow(mgs),ncol=3))
  for (step in 1:nrow(mgs)){
    for (i in 1:2) {
      if(mgs[step,i]<0){
        mgs[step,(i+3)]<- as.numeric(fac[hc$order][(-1)*mgs[step,i]])
      } else {
        mgs[step,(i+3)]<- mgs[mgs[step,i],6]
      }
    }
    if (mgs[step,4]==mgs[step,5]) { mgs[step,6]<- mgs[step,5] } else { mgs[step,6] <-8 }
  }
  
  #an additional matrix is created (dndr_gram)
  #dndr_gram includes info concerning x and y position for each merge
  dndr_gram<-matrix(NA, ncol=3, nrow=nrow(mgs))
  colnames(dndr_gram)<- c("x0","x1","height")
  dndr_gram[,3] <- hc$height
  
  for(step in 1:nrow(mgs)){
    if(mgs[step,1]<0 & mgs[step,2]<0) {
      dndr_gram[step,1]<-mgs[step,1]* (-1)
      dndr_gram[step,2]<-mgs[step,2]* (-1)
    } else if (mgs[step,1] * mgs[step,2] < 0) {
      min <- which(mgs[step,1:2]<0)
      dndr_gram[step,min]<-mgs[step,min]* (-1)
      plu <- which(mgs[step,1:2]>0)
      dndr_gram[step,plu] <- mgs[mgs[step,plu],3] 
    } else {
      dndr_gram[step,1] <- mgs[mgs[step,1],3]
      dndr_gram[step,2] <- mgs[mgs[step,2],3]
    }
  }
  
  #---draw plot---------------------------#
  
  #define y limits & length of leaves
  if (hang <0) { hang = 0 }
  lv_len <- hang * (max(hc$height)-min(hc$height))
  yrange <- 0.1*(max(hc$height)-min(hc$height))
  
  auto.ax.min <- min(hc$height)-lv_len-(0.8*yrange*window.adj)
  auto.ax.max <- max(hc$height+(0.2*yrange*window.adj))
  if(is.numeric(ax.low)) { ax.min <- ax.low } else { ax.min <- auto.ax.min }
  if(is.numeric(ax.top)) { ax.max <- ax.top } else { ax.max <- auto.ax.max }
  
  #conditional checks for setting plot & axis limits
  if(ax.min >= ax.max) {
    ax.min <- auto.ax.min
    ax.max<- auto.ax.max
  }
  
  #plot a dendrogram empty canvas. Legend will be also added.
  plot(1:(nrow(mgs)+3),1:(nrow(mgs)+3),xlab="", ylab="Height", type="n", main=main,
       ylim=c(ax.min,ax.max), axes = F)
  
  ax.low.scale <- ((round(min(hc$height)*0.1, digits=0))*10)-10
  ax.hi.scale <- ((round(max(hc$height)*0.1, digits=0))*10)+10
  if (is.numeric(ax.low) & ax.low < ax.top) {
    ax.low.scale <- round(ax.low, digits=0)
  }
  ax <- seq(ax.low.scale,ax.hi.scale,by=tick.at)
  axis(2, at=ax)
  
  groups <- factor(levels(fac),levels=levels(fac))
  legend("topright",as.vector(groups), pch=15, col=(as.numeric(groups)), bty="n", cex=lab.cex)
  
  #draw vertical segments and labels
  for (step in 1:nrow(mgs)) {
    for (i in 1:2) {
      x<- mgs[step,i]
      if(x<0) {
        segments(x*(-1),dndr_gram[step,3],x*(-1),dndr_gram[step,3]-lv_len,
                 col = if (color == F) { "black" } else { mgs[step,6]}, lwd=lwd)
        text((-1)*x,dndr_gram[step,3]-min((1.4*lv_len),50),hc$labels[hc$order][x*(-1)],
             adj=c(1,0.5), srt=90, cex=lab.cex, col=as.numeric(fac[hc$order][x*(-1)]))
      } else {
        x12<- mgs[x,3]
        y1 <- dndr_gram[step, 3]
        y2 <- dndr_gram[x,3]
        if (color != T) { colr<-"black" } else { colr<-mgs[step,6]}
        segments(x12,y1,x12,y2, col=colr, lwd=lwd)
      }
    }
  }
  #draw horizontal segments
  for (step in 1:nrow(dndr_gram)) {
    if (color != T) { colr<-"black" } else { colr<-mgs[step,6] }
    segments(dndr_gram[step,1],dndr_gram[step,3],dndr_gram[step,2],dndr_gram[step,3],col=colr,lwd=lwd)
  }
}
