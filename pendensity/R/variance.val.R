variance.val <- function(base.den,var.par,weight,K,x,list.len,Z,x.factor,y.val=NULL) {
  M <- floor(K/2)+1
  if(is.null(x)) {
    C <- (diag(c(weight))-tcrossprod(c(weight)))[,-M]
    varr <- t(base.den)%*%C%*%var.par%*%t(C)%*%base.den
    sd <- sqrt(diag(varr))
    return(sd)
  }
  else {
    help.env <- new.env()
    for(i in 1:list.len) {
      ZZ1 <- kronecker(diag(1,K-1),x.factor[i,])
      name2 <- paste("base",i,sep="")
      C.bold <- (diag(weight[i,])-tcrossprod(weight[i,]))[,-M]
      if(!is.null(y.val))
        {
          name3 <- paste("base.y.val",i,sep="")
          base <- get(name3,envir=base.den)
          #C.bold <- (diag(weight[i,])-tcrossprod(weight[i,]))[,-M]
          varr <- t(base)%*%C.bold%*%t(ZZ1)%*%var.par%*%ZZ1%*%t(C.bold)%*%base
          sd <- sqrt(diag(varr))
          assign(paste("sd.cal.y.val",i,sep=""),sd,help.env)
        }
      base <- get(name2,envir=base.den)
      #C.bold <- (diag(weight[i,])-tcrossprod(weight[i,]))[,-M]
      varr <- t(base)%*%C.bold%*%t(ZZ1)%*%var.par%*%ZZ1%*%t(C.bold)%*%base
      sd <- sqrt(diag(varr))
      assign(paste("sd.cal",i,sep=""),sd,help.env)
    }
    return(help.env)
  }
}
