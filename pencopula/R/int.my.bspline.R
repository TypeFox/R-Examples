 int.my.bspline <- function(help.env) {

  stand.num <- c()
  len.k <- length(get("knots.val",help.env)$val)
  if(is.vector(get("base.den",help.env))) base.den <- matrix(base.den,1,length(get("base.den",help.env)))
 
  len.b <- dim(get("base.den",help.env))[2]
  q <- get("q",help.env)#-1 #q ist als order hinterlegt, brauche hier den grad!
  
  knots.val <- get("knots.val",help.env)
      
  #piecewise polynomial calculation

  len.k <- length(knots.val$val)   
      
  #generate help-sequences for calculation
      
  y.all.help <- c()
  for(j in 1:(len.k-1)) {
    help.seq <-  seq(knots.val$val[j],knots.val$val[j+1],length=(q))
    assign(paste("y.help",j,sep=""),help.seq,envir=help.env)
    y.all.help <- c(y.all.help,help.seq)
  }
  
  y.all.help <- unique(y.all.help)

  base.help <- bsplineS(y.all.help,breaks=knots.val$val,norder=get("q",help.env))
  
  for(j in 1:(len.k-1)) {
    list <- which(get("y",help.env)>=knots.val$val[j] & get("y",help.env)<=knots.val$val[j+1])
    assign(paste("y.list",j,sep=""),list,envir=help.env)
    assign(paste("y.part",j,sep=""),get("y",help.env)[list],envir=help.env)
    for(i in 1:(dim(get("base.den",help.env))[2])) { 
      assign(paste("base.part",j,i,sep=""),get("base.den",help.env)[list,i],envir=help.env)
    }
  }  

  #for (i in 1:(len.k-(q-1))) {
  for(i in 1:(len.k-1)) {
    compare <- get(paste("y.help",i,sep=""),envir=help.env)
    list <- which(y.all.help%in%compare)
    for(j in 1:(dim(base.help)[2])) {
      assign(paste("y.base.help",i,j,sep=""),base.help[list,j],envir=help.env)
      assign(paste("y.list.help",i,j,sep=""),list,envir=help.env)
    }
  }
  
  #search the relevant points for calculations und calculate the polynomial-coefficients
  
  q <- q-1 

  for(i in 1:(len.k-1)) {
    y.vec <- c()
    for(j in 1:(dim(base.help)[2])) {
 
      if(q>=0) y.vec <- c(knots.val$val[i])
      if(q>=1) y.vec <- c(y.vec,knots.val$val[i+1])
      if(q>=2) y.vec <- seq(y.vec[1],y.vec[2],length=3)
      if(q>=3) y.vec <- seq(y.vec[1],y.vec[3],length=4)
      if(q>=4) y.vec <- seq(y.vec[1],y.vec[4],length=5)

      assign(paste("y.vec",i,sep=""),y.vec,envir=help.env)
 
      assign(paste("coef",i,".",j,sep=""),(solve(outer(y.vec,0:q,"^"))%*%(get(paste("y.base.help",i,j,sep=""),envir=help.env))),envir=help.env)

    }
  }
      #calculate the integrals and coefficients for standardisation of the splines at the borders
  INT <- matrix(0,dim(base.help)[2],len.k-1)
    
  for(i in 1:(len.k-1)) {

    for(j in 1:(dim(base.help)[2])) {

      y2 <- knots.val$val[i+1]
      y1 <- knots.val$val[i]
      coef <- get(paste("coef",i,".",j,sep=""),envir=help.env)

      y2 <- 1/(1:(q+1))*y2^(1:(q+1))
      y1 <- 1/(1:(q+1))*y1^(1:(q+1))

      INT[j,i] <- sum(coef*y2)-sum(coef*y1)
    }
  }
  assign("INT",INT,help.env)
  if(q==2) {
    list.a <- c(1,dim(INT)[1])
    list.b <- c(1,dim(INT)[2])
    INT.help <- INT[-list.a,]
    INT.help <- INT.help[,-list.b]
    if(is.vector(INT.help)) INT.help <- (1/INT.help) else INT.help <- 1/rowSums(INT.help)
  }
  else  INT.help <- 1/rowSums(INT)
 
  assign("stand.num",INT.help,help.env)
}
