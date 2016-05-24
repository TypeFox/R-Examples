#require(graphics)
## a smallish simple dendrogram

colorDendClass<-function(dend, class)
{
  gI<-0
  mycols=class
  colLab <- function(n,mycols) {
    
    #a <- attributes(n)
    #attr(n, "midpoint") <-c(a$midpoint, 5)
    
    if(is.leaf(n)) {
      a <- attributes(n)
      gI<<-gI+1
      attr(n, "nodePar") <-
        c(a$nodePar, list(lab.col = mycols[gI], lab.cex = 0.5,pch=20,col=mycols[gI]))
      attr(n, "edgePar") <-
        c(a$edgePar, list(col = mycols[gI], lwd=1))
      
      
    }
    return(n)
  }
    
  #op=par();
  dL <- dendrapply(dend, function(n) colLab(n,mycols))
  
  plot(dL,horiz=F,center=FALSE) ## --> colored labels!
  #par(op)
}
