negent3D <-
function(y,nstart=1,m=100,...) {
     f <- function(thetas) {
          cs <- cos(thetas)
          sn <- sin(thetas)
          negent(y%*%c(cs[1],sn[1]*c(cs[2],sn[2])))
     }
     tt <- NULL
     nn <- NULL
     for(i in 1:nstart) {
          thetas <- runif(3)*pi
          o <- optim(thetas,f,method='SANN',control=list(fnscale=-1),...)     
          tt <- rbind(tt,o$par)
          nn <- c(nn,o$value)
     }		
     i <- imax(nn)  # The index of best negentropy
     cs<-cos(tt[i,])
     sn<-sin(tt[i,])
     g.opt <- c(cs[1],sn[1]*cs[2],sn[1]*sn[2])
     g.opt <- cbind(g.opt,c(-sn[1],cs[1]*cs[2],sn[2]*cs[1]))
     g.opt <- cbind(g.opt,c(0,-sn[2],cs[2]))
     x <- y%*%g.opt[,2:3]
     n2 <- negent2D(x,m=m)
     g.opt[,2:3] <- g.opt[,2:3]%*%n2$vectors
     list(vectors=g.opt,values = c(nn[i],n2$values))
}
