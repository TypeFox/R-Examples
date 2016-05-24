smooth12 <- function(x,y,xout,knum=16,std=TRUE){

n = length(y)
nk = ncol(as.matrix(x))
if (nk==1){
  target <- as.matrix(xout)
  otarget1 <- target
  target <- sort(unique(target))
  otarget2 <- target
 
  o <- order(x)
  x <- x[o]
  y <- y[o]
  y <- aggregate(y,by=list(x),"mean")[,2]
  sampvar <- !duplicated(x)
  x <- x[sampvar]

  ylag2 <- c(NA,NA,y[1:(n-2)])
  ylag1 <- c(NA,y[1:(n-1)])
  y1 <- c(y[2:n],NA)
  y2 <- c(y[3:n],NA,NA)
  xlag2 <- c(NA,NA,x[1:(n-2)])
  xlag1 <- c(NA,x[1:(n-1)])
  x1 <- c(x[2:n],NA)
  x2 <- c(x[3:n],NA,NA)
  m1 <- (ylag1-ylag2)/(xlag1-xlag2)
  m2 <- (y-ylag1)/(x-xlag1)
  m3 <- (y1-y)/(x1-x)
  m4 <- (y2-y1)/(x2-x1)
  slope <- (abs(m4-m3)*m2 + abs(m2-m1)*m3)/(abs(m4-m3) + abs(m2-m1))
  slope[is.na(m1)&is.na(m2)] <- m3[is.na(m1)&is.na(m2)]
  slope[is.na(m3)&is.na(m4)] <- m2[is.na(m3)&is.na(m4)]
  sampvar <- is.na(m1)&!is.na(m2)
  slope[sampvar] <-  (abs(m4[sampvar]-m3[sampvar])*m2[sampvar] + m3[sampvar])/(abs(m4[sampvar]-m3[sampvar]) + 1)
  sampvar <- is.na(m4)&!is.na(m3)
  slope[sampvar] <- (m2[sampvar] + abs(m2[sampvar]-m1[sampvar])*m3[sampvar])/(1 + abs(m2[sampvar]-m1[sampvar]))
  mvar <- c(m1,m2,m3,m4)
  mvar <- mvar[!is.na(mvar)]
  slope <- ifelse(is.na(slope), mean(mvar), slope)
  xdata <- data.frame(x,y,slope)

# merge and sort target and x data so that x's are between target points
  nt = length(target)
  z <- array(0,dim=nt)
  tdata <- data.frame(target)
  names(tdata) <- "x"
  xdata <- merge(xdata,tdata,by="x",all=TRUE)
  x1 <- xdata$x[!is.na(xdata$y)]
  t1 <- xdata$slope[!is.na(xdata$y)]
  y1 <- xdata$y[!is.na(xdata$y)]
  xsum <- cumsum(!is.na(xdata$y))
  xdata$x1 <- ifelse(xsum>0&xsum<n, x1[xsum], NA)
  xdata$t1 <- ifelse(xsum>0&xsum<n, t1[xsum], NA)
  xdata$y1 <- ifelse(xsum>0&xsum<n, y1[xsum], NA)
  xdata$x2 <- ifelse(xsum>0&xsum<n, x1[xsum+1], NA)
  xdata$t2 <- ifelse(xsum>0&xsum<n, t1[xsum+1], NA)
  xdata$y2 <- ifelse(xsum>0&xsum<n, y1[xsum+1], NA)
  xdata$target <- xdata$x


#loader 
  lam <- (xdata$x - xdata$x1)/(xdata$x2-xdata$x1)
  phi1 <- ((1-lam)^2)*(1 + 2*lam)
  phi2 <- (lam^2)*(3-2*lam)
  psi1 <- lam*((1-lam)^2)
  psi2 <- -(lam^2)*(1-lam)
  ytarget <- phi1*xdata$y1 + phi2*xdata$y2 + (xdata$x2-xdata$x1)*(psi1*xdata$t1 + psi2*xdata$t2)
  ytarget <- ytarget[xdata$target%in%target]
  if (target[1]==min(x)){ytarget[1] = y[1]}
  if (target[nt]==max(x)){ytarget[nt] = y[n]}

  data1 <- data.frame(otarget1)
  names(data1) <- "target"
  data2 <- data.frame(otarget2,ytarget)
  names(data2) <- c("target","ytarget")

  data1 <- merge(data1,data2,by="target",sort=FALSE)
  o1 <- order(xout)
  o2 <- order(data1$target)
  data1 <- data1[o2,]
  data1[o1,] <- data1

  names(data1) <- c("x","y")
  yout <- data1$y
}
if (nk==2){
  if (std==TRUE){
    vx <- diag(1/apply(x,2,sd))
    x <- data.frame(as.matrix(x)%*%vx)
    xout <- data.frame(as.matrix(xout)%*%vx)
  }

  fit <- nn2(x,xout,k=knum)
  maxd = apply(fit$nn.dists,1,max)
  yout <- array(0,dim=nrow(xout))
  wsum <- yout
  for (j in seq(1,knum)){
    w <- ifelse(fit$nn.dists[,j]>0.000000001,((maxd - fit$nn.dists[,j])^2)/(maxd*fit$nn.dists[,j]),1)
    yout <- yout + w*y[fit$nn.idx[,j]]
    wsum <- wsum+w
  }
  yout <- yout/wsum
}

  return(yout)
}

 
