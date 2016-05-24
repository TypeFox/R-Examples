Derv2 <- function(penden.env,temp=FALSE,lambda=NULL) {
  correct <- 0.00001 * diag(get("DD",penden.env))
  if(!temp) {
    Fy <- kronecker(get("tilde.PSI.d.D",penden.env) %*% get("ck.val",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))
    assign("Derv2.pen",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)-get("lambda",penden.env)*(get("DDD.sum",penden.env))-correct),penden.env)
    assign("Derv2.cal",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)),penden.env)
  }
  if(temp) {
    if(is.null(lambda)) lambda <- get("lambda",penden.env)
    Fy <- kronecker(get("tilde.PSI.d.D",penden.env) %*% get("ck.val.temp",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))
    assign("Derv2.cal.temp",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)),penden.env)
    assign("Derv2.pen.temp",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)-lambda*(get("DDD.sum",penden.env))-correct),penden.env)
  }
}
