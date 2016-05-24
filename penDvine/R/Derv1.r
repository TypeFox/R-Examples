Derv1 <- function(penden.env,temp=FALSE,lambda=NULL) {
  if(!temp) {
    assign("Derv1.pen",matrix(colSums(get("tilde.PSI.d.D",penden.env)/kronecker(get("f.hat.val",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))),get("DD",penden.env),1)-get("lambda",penden.env)*get("DDD.sum",penden.env)%*%get("ck.val",penden.env),penden.env)
  }
  else {
    if(is.null(lambda)) lambda <- get("lambda",penden.env)
    assign("Derv1.pen.temp",matrix(colSums(get("tilde.PSI.d.D",penden.env)/kronecker(get("f.hat.val.temp",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))),get("DD",penden.env),1)-lambda*get("DDD.sum",penden.env)%*%get("ck.val.temp",penden.env),penden.env)
  } 
}
