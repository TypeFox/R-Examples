 my.bspline <- function(h,q,knots,y,K,plot.bsp,typ) {
  len.k <- length(knots)
  
  base.den <- bsplineS(y,breaks=knots,norder=q)
  len.b <- dim(base.den)[2]

  knots.val <- list()
  knots.val$val <- knots
  
  #Normierung

  if(typ==1) {
    stand.num <- c()
    min.knot <- min(knots)
    max.knot <- max(knots)
    knots.val$help <- unique(c(seq(min.knot-(q-1)*h,min.knot-h,by=h),knots,seq(max.knot+h,max.knot+(q-1)*h,by=h)))
    for(i in 1:(length(knots.val$help-q))) stand.num[i] <- (q/(knots.val$help[i+q]-knots.val$help[i]))
    INT <- NULL
  }

  #integration
  
  if(typ==2 | typ==3) {
    help.env <- new.env()
    assign("base.den",base.den,help.env)
    assign("knots.val",knots.val,help.env)
    assign("y",y,help.env)
    assign("q",q,help.env)
    int.my.bspline(help.env)
    stand.num <- get("stand.num",help.env)
    INT <- get("INT",help.env)
  }
  
  if(typ==1|typ==2|typ==3) {
    if(q==3) {
      base.den <- base.den[,-c(1,dim(base.den)[2])]
      len.b <- len.b -2
    }
    for(j in 1:len.b) base.den[,j] <- base.den[,j]*stand.num[j]
  }

  #printing
  if(plot.bsp==TRUE) {
    x1 <- c()
    fx1 <- c()
    fac1 <- c()
    for(i in 1:len.b) {
      help <- y
      x1 <- c(x1,help)
      fx1 <- c(fx1,base.den[,i])
      if(i<=9) fac1 <- c(fac1,rep(paste("Bspline-No.",0,i,sep=""),length(help)))
      if(i>9)fac1 <- c(fac1,rep(paste("Bspline-No.",i,sep=""),length(help)))
    }
    datafr1 <- data.frame(x1,fx1,fac1)
    help2 <- xyplot(fx1~x1,groups=fac1,type=c("p"),scales=list(x='free',y='free'),ylab="",as.table=TRUE,data=datafr1,main=paste("Distribution of ",len.b," B-Splines, degree ",q-1),auto.key=list(space="bottom",columns=6),xlab="y")
    print(help2)
  }
  return(list(base.den=base.den,stand.num=stand.num,knots.val=knots.val,K=K,INT=INT))
}

