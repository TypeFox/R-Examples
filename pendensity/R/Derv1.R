#first derivative - with and without covariables
Derv1 <- function(penden.env){
  ck.temp <- get("ck.temp",penden.env)
  M <- get("M",penden.env)
  K <- get("K",penden.env)
  N <- get("N",penden.env)
  base <- get("base",penden.env)
  Z.index <- get("Z.index",penden.env)
  len.x.fac <- get("len.x.fac",penden.env)

  for(i in 1:len.x.fac) {
    if(is.null(get("x",penden.env))) assign("C.bold",t((diag(c(ck.temp))-crossprod(ck.temp))[,-M]),penden.env)
    else assign(paste("C.bold",i,sep=""),t((diag((ck.temp[i,]))-tcrossprod(ck.temp[i,]))[,-M]),penden.env)
  }
  index <- matrix(1:get("n",penden.env))

  #penalty
  #beta.val.help <- c(get("beta.val",penden.env)[1:(N*(M-1))],get("beta.val",penden.env)[(N*M+1):(K*N)])
  #penalty <- get("lambda0",penden.env)*get("Dm",penden.env)%*%c(get("beta.val",penden.env)[1:(N*(M-1))],get("beta.val",penden.env)[(N*M+1):(K*N)])

  #built the first derivative for each y_i, depending on the base choosen
  if(base=="bspline") {
    if(get("x.null",penden.env)) Derv1.cal <- apply(index,1,function(i,base.den,f.hat.val) get("C.bold",penden.env)%*%base.den[,i]/f.hat.val[i] ,get("base.den",penden.env),get("f.hat.val",penden.env))
    else  Derv1.cal <- apply(index,1,function(i,base.den,f.hat.val,Z,Z.index,x.factor) outer(x.factor[Z.index[i],],get(paste("C.bold",Z.index[i],sep=""),penden.env)%*%(base.den[,i]/f.hat.val[i]),"*"),get("base.den",penden.env),get("f.hat.val",penden.env),get("Z",penden.env),Z.index,get("x.factor",penden.env))
  }  
  if(base=="gaussian") {
    if(get("x.null",penden.env)) Derv1.cal <- apply(index,1,function(i,base.den,f.hat.val,Stand.abw) get("C.bold",penden.env)%*%diag(1/Stand.abw)%*%base.den[,i]/f.hat.val[i] ,get("base.den",penden.env),get("f.hat.val",penden.env),get("Stand.abw",penden.env))
    else Derv1.cal <- apply(index,1,function(i,base.den,f.hat.val,Stand.abw,Z,Z.index,x.factor) outer(x.factor[Z.index[i],],get(paste("C.bold",Z.index[i],sep=""),penden.env)%*%diag(1/Stand.abw)%*%(get("base.den",penden.env)[,i]/f.hat.val[i]),"*"),get("base.den",penden.env),get("f.hat.val",penden.env),get("Stand.abw",penden.env),get("Z",penden.env),Z.index,get("x.factor",penden.env))
  }
  #Derv1.pen <- rowSums(Derv1.cal)-penalty
  return(list(Derv1.pen=rowSums(Derv1.cal)-get("lambda0",penden.env)*get("Dm",penden.env)%*%c(get("beta.val",penden.env)[1:(N*(M-1))],get("beta.val",penden.env)[(N*M+1):(K*N)]),Derv1.cal=Derv1.cal,f.hat.val=get("f.hat.val",penden.env)))  
}
 
