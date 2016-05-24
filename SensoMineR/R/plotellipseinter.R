"plotellipseinter" <- function(mat,alpha=0.05,coord=c(1,2),nbgroup=1,moy=TRUE,eig,cex=1,color=NULL,title=NULL){

#################################################################
"ellipse2" <- function(loc, cov,alpha)
      {
            A <- cov
            detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
            dist <- sqrt(qchisq(1-alpha, 2))
            ylimit <- sqrt(A[2, 2]) * dist
            y <- seq( - ylimit, ylimit, 0.01 * ylimit)
            sqrt.discr <- sqrt(detA/A[2, 2]^2 * abs(A[2, 2] * dist^2 - y^2))
            sqrt.discr[c(1, length(sqrt.discr))] <- 0
            b <- loc[1] + A[1, 2]/A[2, 2] * y
            x1 <- b - sqrt.discr
            x2 <- b + sqrt.discr
            y <- loc[2] + y
            return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
      }
#################################################################

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
if (moy ==TRUE) {
  matP=cbind.data.frame(mat$moy$P[,coord],mat$moy$P[,ncol(mat$moy$P)])
  matPJ=cbind.data.frame(mat$moy$PJ[,coord],mat$moy$PJ[,ncol(mat$moy$PJ)])
  matsimul=cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
}
if (moy == FALSE) {
  matmoyP=cbind.data.frame(mat$moy$P[,coord],mat$moy$P[,ncol(mat$moy$P)])
  matmoyPJ=cbind.data.frame(mat$moy$PJ[,coord],mat$moy$PJ[,ncol(mat$moy$PJ)])
  matmoysimul=cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
  matP=cbind.data.frame(mat$partiel$P[,coord],mat$partiel$P[,ncol(mat$partiel$P)])
  matPJ=cbind.data.frame(mat$partiel$PJ[,coord],mat$partiel$PJ[,ncol(mat$partiel$PJ)])
  matsimul=cbind.data.frame(mat$partiel$simul[,coord],mat$partiel$simul[,ncol(mat$partiel$simul)])
}  
nbp <- nrow(matP)
nbprod <- nbp/nbgroup
coord.ellipse.a.tracer <- matrix(0,402,2*nbp)

p <- 2
nbjuge <-  nrow(matPJ)/nrow(matP)
nbsimul <-  nrow(matsimul)/nrow(matP)

for (i in 1:nbp){
  VX <- var(matsimul[((i-1)*nbsimul+1):(i*nbsimul),1:2])
  coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse2(t(matP[i,1:2]),VX,alpha)
}

minx <- min(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
maxx <- max(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
miny <- min(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)
maxy <- max(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)

  plot(0, 0, xlab = paste("Dim ",coord[1]," (",eig[coord[1],2],"%)",sep=""), ylab = paste("Dim ",coord[2]," (",eig[coord[2],2],"%)",sep=""), xlim = c(minx*1.05,maxx*1.05), ylim = c(1.05*miny,1.05*maxy), col = "white", asp=1)
  if (is.null(title)){
    if (moy==TRUE) title(main = "Confidence ellipses for the mean points")
    if (moy==FALSE) title(main = "Confidence ellipses for the partial points")
  } else {title(main=title)}
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  if (moy==FALSE){
    points(matmoyP[,1], matmoyP[,2],cex=0.8*cex,col=color[1:nbprod],pch=15)
    text( matmoyP[,1], matmoyP[,2], matmoyP[,ncol(matmoyP)], cex = 0.8*cex, pos = 4, offset = 0.2,col=color[1:nbprod])
  }
  if (moy==TRUE) text( matP[,1], matP[,2], matP[,ncol(matP)], cex = 0.8*cex, pos = 4, offset = 0.2,col=color[1:nbprod])
  for (j in 1:nbgroup){
    for (i in 1:nbprod) {
      points(matP[(j-1)*nbprod+i,1], matP[(j-1)*nbprod+i,2],cex=0.8*cex,col=color[i],pch=20)
      if (moy==FALSE) lines(c(matP[(j-1)*nbprod+i,1],matmoyP[i,1]),c( matP[(j-1)*nbprod+i,2],matmoyP[i,2]),col=color[i],lty=j)
      lines(coord.ellipse.a.tracer[,(1+2*((i+(j-1)*nbprod)-1)):(2*(i+(j-1)*nbprod))],col=color[i],lty=j)
    }
  }
}
