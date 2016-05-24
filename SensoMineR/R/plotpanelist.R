"plotpanelist" <- function(mat,coord=c(1,2),name=FALSE,eig,cex=1,color=NULL)   {

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
mat=cbind.data.frame(mat[,coord],mat[,(ncol(mat)-1):ncol(mat)])
nbprod <- length(levels(mat[,ncol(mat)-1]))
lab <- mat[1:nbprod,ncol(mat)-1]
nbjuge <-  length(levels(mat[,ncol(mat)]))-1
nbpoint=nbprod*(nbjuge+1)

minx <- min(mat[1:nbpoint,1],na.rm=TRUE)
maxx <- max(mat[1:nbpoint,1],na.rm=TRUE)
miny <- min(mat[1:nbpoint,2],na.rm=TRUE)
maxy <- max(mat[1:nbpoint,2],na.rm=TRUE)

  plot(0, 0, xlab = paste("Dim ",coord[1]," (",eig[coord[1],2],"%)",sep=""), ylab = paste("Dim ",coord[2]," (",eig[coord[2],2],"%)",sep=""), xlim = c(minx*1.05,maxx*1.05), ylim = c(1.05*miny,1.05*maxy), col = "white", asp=1)
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  points(mat[1:nbprod,1],mat[1:nbprod,2],cex=cex*1.2,col=color[1:nbprod],pch=15)
  text( mat[1:nbprod,1], mat[1:nbprod,2], labels=lab, cex = cex*0.8, pos = 4, offset = 0.2,col=color[1:nbprod])
  if (name==FALSE) for (i in 1:nbjuge) points(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],cex=cex*0.8,col=color[1:nbprod],pch=20)
  if (name==TRUE) for (i in 1:nbjuge) text(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],labels=mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),ncol(mat)],cex=cex*0.6,col=color[1:nbprod],pch=20)
  title(main = "Individual description")

}
