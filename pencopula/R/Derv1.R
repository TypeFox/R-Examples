 Derv1 <- function(penden.env) {
  assign("Derv1.pen",matrix(colSums(get("tilde.PSI.d.D",penden.env)/kronecker(get("tilde.PSI.d.D",penden.env)%*%get("ck.val",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))),get("DD",penden.env),1)-get("DDD.sum",penden.env)%*%get("ck.val",penden.env),penden.env)
}
