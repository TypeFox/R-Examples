HIndivabs <- function (data, catch="Food", hand="Hand", indiv="Indiv", RightHand = "R", LeftHand = "L"
                       , col = 1:length(levels(data[[catch]])), ylab = "Absolute handedness index"
                       , main="Hand preference regarding to the performed task by each individual", cex.main=1
                       , legend.text = FALSE, beside = TRUE, ylim = c(0,1), vlines = TRUE, hlines = TRUE
                       , legendlocation=TRUE, cex=1, pt.cex=2, pch=15, savetable = FALSE, file = "HIperIndivabs.csv")
{
  for (i in 1:nlevels(data[[catch]])) {
      seldata<- data[data[[catch]]==levels(data[[catch]])[i],]
      Tab<- table(seldata[[indiv]], seldata[[hand]])
      NewTab<-as.data.frame.matrix(Tab)
      ifelse (is.null(NewTab[[RightHand]]) == TRUE, HITab<-(-NewTab[[LeftHand]])/NewTab[[LeftHand]], ifelse (is.null(NewTab[[LeftHand]]) == TRUE, HITab<-NewTab[[RightHand]]/NewTab[[RightHand]], HITab<-(NewTab[[RightHand]]-NewTab[[LeftHand]])/(NewTab[[RightHand]]+NewTab[[LeftHand]]))) #Handedness index            
      HITab[which(HITab<0)]<-HITab[which(HITab<0)]*-1
      if("HIabsperIndiv" %in% ls() == FALSE) {HIabsperIndiv<-c()} else {}
      HIabsperIndiv<-cbind(HIabsperIndiv,HITab)
  }
  HIabsperIndiv<-t(HIabsperIndiv)
  colnames(HIabsperIndiv)<-levels(data[[indiv]])
  rownames(HIabsperIndiv)<-levels(data[[catch]])
  
  graph<-as.matrix(HIabsperIndiv)
  barplot(graph, beside = beside, ylab=ylab, main=main, cex.main=cex.main, legend.text = legend.text, col=col, ylim=ylim)

  # Vertical lines plot
  if (nlevels(data[[indiv]])>1) {
      if (vlines == TRUE) {
        abline(v=(seq(nlevels(data[[catch]])+1.5,nlevels(data[[indiv]])*(nlevels(data[[catch]])+1),by=nlevels(data[[catch]])+1)), lty=3)
      } else {
        }
  } else {
    }
  
  # Horizontal lines plot
  if (hlines == TRUE) {
      ya<-rep(0,nlevels(data[[indiv]]))
      yb<-rep(0,nlevels(data[[indiv]]))
      y<-cbind(ya,yb)
      xa<-seq(1,(nlevels(data[[indiv]])-1)*(nlevels(data[[catch]])+1)+1,by=nlevels(data[[catch]])+1)
      xb<-seq(1+nlevels(data[[catch]]),nlevels(data[[indiv]])*(nlevels(data[[catch]])+1),by=nlevels(data[[catch]])+1)
      x<-cbind(xa,xb)
      for (i in 1:nlevels(data[[indiv]])) {
        lines (x[i,],y[i,])
      }
  } else {
    }

  #Legend
  if (legendlocation == TRUE) {
      message("Click where you want to place the legend")
      legendplace <- locator(1)
      legend(legendplace$x,legendplace$y,as.vector(levels(data[[catch]])),col=col,bty="n",pch=pch, cex=cex, pt.cex=pt.cex)
  } else {
    }

  if (savetable == "csv") {write.csv(HIabsperIndiv, file = file)} else {}
  if (savetable == "csv2") {write.csv2(HIabsperIndiv, file = file)} else {}    
  HIabsperIndiv
}

