my.bspline <- function(h,q,knots.val,y,K,plot.bsp) {
  help.degree <- 0
  if(q>2) help.degree <- q-2
 
  knots.spline <- knots.val$val
  len.k <- length(knots.spline)
  knots.spline <- c(knots.spline[1]-6*h,knots.spline[1]-5*h,knots.spline[1]-4*h,knots.spline[1]-3*h,knots.spline[1]-2*h,knots.spline[1]-h,knots.spline,knots.spline[len.k]+h,knots.spline[len.k]+2*h,knots.spline[len.k]+3*h,knots.spline[len.k]+4*h,knots.spline[len.k]+5*h,knots.spline[len.k]+6*h)

  knots.add.l <- len.k+2*(q-1)
  list.knots.add <- (1:knots.add.l)+(6-(q-1))
  knots.val$help <- knots.spline[list.knots.add]
  knots.val <- list(dis=knots.val$dis,val=knots.val$val,all=knots.spline,help=knots.val$help)
  
  base.den <- t(bsplineS(y,breaks=knots.spline,norder=q))
  base.den2 <- t(bsplineS(y,breaks=knots.spline,norder=(q+1)))

  K <- K+help.degree
  K1 <- K+1
  base.den <- base.den[7:(K+6),]
  base.den2 <- base.den2[7:(K1+6),]
 
  #Normierung
  
  stand.num <- c()
  for(i in 1:(length(knots.val$help)-q)) {
    stand.num[i] <- (q/(knots.val$help[i+q]-knots.val$help[i]))
  }
  base.den <- stand.num*base.den
  
  len.b <- length(base.den[,1])

  #printing
  if(plot.bsp==TRUE) {
    x1 <- c()
    fx1 <- c()
    fac1 <- c()
    for(i in 1:len.b) {
      help <- y
      x1 <- c(x1,help)
      fx1 <- c(fx1,base.den[i,])
      if(i<=9) fac1 <- c(fac1,rep(paste("Bspline-No.",0,i,sep=""),length(help)))
      if(i>9)fac1 <- c(fac1,rep(paste("Bspline-No.",i,sep=""),length(help)))
    }
    datafr1 <- data.frame(x1,fx1,fac1)
    help2 <- xyplot(fx1~x1,groups=fac1,type=c("p"),scales=list(x='free',y='free'),ylab="",as.table=TRUE,data=datafr1,main=paste("Distribution of ",len.b," B-Splines, degree ",q-1),auto.key=list(space="bottom",columns=6),xlab="y")
    print(help2,more=TRUE)
  } 
  return(list(base.den=base.den,base.den2=base.den2,stand.num=stand.num,knots.val=knots.val,K=K,help.degree=help.degree))
}
