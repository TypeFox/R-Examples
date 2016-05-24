distr.func.help <- function(base,knots,penden.env,q,y) {
  x.help <- c()
  n <- length(y)
  for(j in 1:(length(knots)-1)) {
    x.help <- c(x.help,(knots[j+1]-knots[j])/2+knots[j])
  }
  knots.help <- sort(c(x.help,knots))
  base.help <- my.bspline(y=knots.help, K = (get("K",penden.env)+get("q",penden.env)-1),q=q)$base.den
  len.b <- dim(base.help)[2]
  if(is.vector(base.help)) {
    base.help <- matrix(base.help,length(base.help),1)
    len.b <- 1
  }
  #knots.help <- get("knots.help",penden.env)
   
  knots.list <- which(knots.help%in%knots)
  knots.l <- knots.help[knots.list]
  help.env <- new.env()
  len.k <- length(knots.l)
  #q <- q+1
  y.all.help <- c()

  #help points between knots
  
  #for(i in 1:(len.k-(q-1))) {
  for(i in 1:(len.k-1)) {
    a <- which(knots.help==knots.l[i])
    if(q==2) b <- which(knots.help==knots.l[i+1])
    if(q==1) b <- a+1
    if(q==1) help.seq <- knots.help[c(a,b)]
    if(q==2) help.seq <- knots.help[c(a,(a+b)/2,b)]
    if(q==3) help.seq <- knots.help[c(a,a+(b-a)/3,a+2*(b-a)/3,b)]
    assign(paste("y.help",i,sep=""),help.seq,envir=help.env)
    y.all.help <- c(y.all.help,help.seq)
  }  

  y.all.help <- unique(y.all.help)
 
  #which help points are for which base part
  
  for(i in 1:(len.k-1)) {
    for (j in 1:(len.b)) {
      compare <- get(paste("y.help",i,sep=""),envir=help.env)#i oder j?
      list <- which(knots.help%in%compare)
      #print(list)
      assign(paste("y.base.help",i,".",j,sep=""),base.help[list,j],envir=help.env)
      #print(base.help[list,j])
      assign(paste("y.list.help",i,".",j,sep=""),list,envir=help.env)
    }
  }
  
#############  
  #search the relevant points for calculations und calcute the polynomial-coefficients of each base part

#q <- q-1

for(j in 1:len.b) {
  #print("j")
  #print(j)
    for(i in 1:(len.k-1)) {
      #print("i")
      #print(i)
      y.vec <- c()
      a <- which(knots.help==knots.l[i])
      b <- which(knots.help==knots.l[i+1])
      if(q>=0) y.vec <- knots.help[a]
      if(q>=1) y.vec <- c(y.vec,knots.help[b])
      if(q==2) y.vec <- c(y.vec[1],knots.help[(which(knots.help==y.vec[1])+which(knots.help==y.vec[2]))/2],y.vec[2])
      if(q==3) y.vec <- knots.help[c(a,a+(b-a)/3,a+2*(b-a)/3,b)]
      if(q>=4) y.vec <- seq(y.vec[1],y.vec[4],length=5)
      assign(paste("y.vec",i,".",j,sep=""),y.vec,envir=help.env)
      assign(paste("coef",i,".",j,sep=""),(solve(outer(y.vec,0:q,"^"))%*%(get(paste("y.base.help",i,".",j,sep=""),envir=help.env))),envir=help.env)
    }
  }

  yi <- NULL

  func.env <- new.env()
  sum <- val <- 0

  knots <- knots.l
  obj<-NULL
  rm(obj)
  for(j in 1:len.b) {
    for(i in 1:(len.k-1)) {
      term <- paste("(",poly.part(i,j,knots,help.env,q,poly=TRUE),")",sep="")
      assign(paste("distr.func",i,".",j,sep=""),paste("obj <-function(x){",term,"}"),envir=func.env)
    }
  }
  
  y.val <- rep(0,n)
  int.base <- matrix(0,length(y),dim(base.help)[2])
  
  for(j in 1:len.b) {
    INT <- 0
    for(i in 1:(len.k-1)) {
      y.val <- rep(0,n)
      funcy <- get(paste("distr.func",i,".",j,sep=""),envir=func.env)
      eval(parse(text=funcy))
      if(i<(len.k-1)) list <- which(knots.l[i] <=y & knots.l[i+1] > y)
      if(i==(len.k-1)) list <- which(knots.l[i] <=y & knots.l[i+1] >= y)
      if(i==1) int.base[list,j] <- obj(y[list]) + int.base[list,j]
      if(i>1) int.base[list,j] <- obj(y[list]) + int.base[list,j] + INT
      INT <- obj(knots.l[i+1]) + INT
    }
  }
  return(int.base)
}

