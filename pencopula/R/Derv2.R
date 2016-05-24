Derv2 <- function(penden.env,temp=FALSE) {
  if(!temp) {
    if(get("base",penden.env)=="Bernstein") correct <- 0.00001 * diag(get("DD",penden.env)) else correct <- 0
    Fy <- kronecker(get("tilde.PSI.d.D",penden.env) %*% get("ck.val",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))
    for(j in 1:get("p",penden.env)) {
      assign(paste("Derv2.pen.",j,sep=""),-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)+(-get("DDD",penden.env)[,,j])-correct,penden.env)
    }
    assign("Derv2.pen",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)+(-get("DDD.sum",penden.env)))-correct,penden.env)
    assign("Derv2.cal",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)),penden.env)
  }
  if(temp) {
    Fy <- kronecker(get("tilde.PSI.d.D",penden.env) %*% get("ck.val.temp",penden.env), matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))
    for(j in 1:get("p",penden.env)) {
      assign(paste("Derv2.pen.temp.",j,sep=""),-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)+(-get("DDD.temp",penden.env)[,,j])-correct,penden.env)
    }
    assign("Derv2.cal.temp",(-crossprod(get("tilde.PSI.d.D",penden.env)/Fy)),penden.env)
  }
} 
