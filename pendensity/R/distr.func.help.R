distr.func.help <- function(obj) {
  y <- obj$values$y
  q <- obj$splines$q
  knots.val <- obj$splines$knots.val
  base.den <- obj$splines$base.den
  base.den2 <- obj$splines$base.den2
  len.b <- length(base.den[,1])
  help.env <- new.env()
  K <- obj$splines$K
  h <- obj$splines$h
  diff.b <- abs(len.b-length(knots.val$help))-1
  help.degree <- obj$splines$help.degree
  
  y.all.help <- c()

  #help points between knots
  
  for(i in 1:(length(knots.val$help)-1)){
    help.seq <-  seq(knots.val$help[i],knots.val$help[i+1],length=(q+1))
    assign(paste("y.help",i,sep=""),help.seq,envir=help.env)
    y.all.help <- c(y.all.help,seq(knots.val$help[i],knots.val$help[i+1],length=(q+1)))
  }  
  #for(i in len.b:(len.b+(q-1))) y.all.help <- c(y.all.help,seq(knots.val$help[i],knots.val$help[i+1],length=(q+1)))

  y.all.help <- unique(y.all.help)
  base.norm <- my.bspline(h,q,knots.val,y=y.all.help,K-help.degree,plot.bsp=FALSE)
  base.help <- base.norm$base.den
  #which help points are for which base part

  for (i in 1:(len.b)) {
    for(j in 1:q) {
      compare <- get(paste("y.help",(i+j-1),sep=""),envir=help.env)
      list <- which(y.all.help%in%compare)
      base.val <- base.help[i,list]
      assign(paste("y.base.help",i,".",j,sep=""),base.val,envir=help.env)
      assign(paste("y.list.help",i,".",j,sep=""),list,envir=help.env)
    }
  }

  #search the relevant points for calculations und calcute the polynomial-coefficients of each base part

  for(i in 1:(len.b)) {
    y.vec <- c()
    for(j in 1:q) {
      if(q>=0) y.vec <- c(knots.val$help[i+j])
      if(q>=1) y.vec <- c(knots.val$help[i+j-1],y.vec)
      if(q>=2) y.vec <- seq(y.vec[1],y.vec[2],length=3)
      if(q>=3) y.vec <- seq(y.vec[1],y.vec[3],length=4)
      if(q>=4) y.vec <- seq(y.vec[1],y.vec[4],length=5)
      assign(paste("y.vec",i,".",j,sep=""),y.vec,envir=help.env)
      assign(paste("coef",i,".",j,sep=""),(solve(outer(y.vec,0:q,"^"))%*%(get(paste("y.base.help",i,".",j,sep=""),envir=help.env))),envir=help.env)
    }
  }
  
  help.env <- cal.int(len.b,q,help.env,knots.val)

  return(help.env)
}
