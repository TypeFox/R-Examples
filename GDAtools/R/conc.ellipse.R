conc.ellipse <- function(resmca,var,sel=1:length(levels(varb)),col=rainbow(length(sel)),axes=c(1,2),cex=0.2) {
  m <- varsup(resmca,var)$coord[,axes]
  m[,1] <- m[,1]*resmca$svd$vs[axes[1]]
  m[,2] <- m[,2]*resmca$svd$vs[axes[2]]
  v <- varsup(resmca,var)$var[1:length(levels(var)),axes]
  classe <- class(resmca)[1] # new
  if(classe=='stMCA') classe=resmca$call$input.mca # new
  if(classe == 'csMCA') { # new
     varb <- var[resmca$call$subcloud]
     wt <-  resmca$call$row.w[resmca$call$subcloud]
     }
  if(classe %in% c('MCA','speMCA','multiMCA')) { # new
     varb <- var
     wt <-  resmca$call$row.w
     }
  c <- vector(length=length(levels(var)))
  for(i in 1:length(c)) {
     temp1 <- matrix(resmca$ind$coord[varb==levels(varb)[i],axes],ncol=2)
     temp1[,1] <- temp1[,1] - m[i,1]
     temp1[,2] <- temp1[,2] - m[i,2]
     temp2 <- wt[varb==levels(varb)[i]]*temp1[,1]*temp1[,2]
     c[i] <- sum(temp2)/sum(wt[varb==levels(varb)[i]])
     }
  g1 <- 0.5*(v[,1]+v[,2])+0.5*sqrt((v[,1]-v[,2])^2+4*c^2)
  g2 <- 0.5*(v[,1]+v[,2])-0.5*sqrt((v[,1]-v[,2])^2+4*c^2)
  sa1 <- 2*sqrt(g1)
  sa2 <- 2*sqrt(g2)
  alph <- atan((g1-v[,1])/c)
  npoints <- 100
  theta <- seq(0, 2 * pi, length=(npoints))
  if(!(is.null(col)) & length(col)!=length(sel))  col <- rainbow(length(sel))
  for(i in 1:length(levels(varb))) {
    if(i %in% sel) {
      colprinc <- col[which(i==sel)]
      x0 <- m[i,1]
      y0 <- m[i,2]
      alpha <- alph[i]
      a <- sa1[i]
      b <- sa2[i]
      x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
      y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
      lines(x, y, type = "l", col=colprinc)
      z1 <- x0 + c(a,b,a,b)*cos(alpha+c(0,pi/2,pi,3*pi/2))
      z2 <- y0 + c(a,b,a,b)*sin(alpha+c(0,pi/2,pi,3*pi/2))
      z <- cbind(z1,z2)
      lightcol=rgb(t(col2rgb(colprinc)),alpha=200,maxColorValue=255)
      lines(z[c(1,3),],col=lightcol,lty=2)
      lines(z[c(2,4),],col=lightcol,lty=2)
      points(x0,y0,pch=19,cex=1,col=colprinc)
      points(resmca$ind$coord[varb==levels(varb)[i],axes],pch=19,cex=cex,col=colprinc)
      text(x0,y0+0.1,levels(varb)[i],col=colprinc)
      }
    }
}
