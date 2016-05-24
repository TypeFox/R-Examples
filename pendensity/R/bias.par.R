bias.par <- function(penden.env) {
  M <- get("M",penden.env)
  N <- get("N",penden.env)
  K <- get("K",penden.env)
  beta.val.help <- c(get("beta.val",penden.env)[1:(N*(M-1))],get("beta.val",penden.env)[(N*M+1):(K*N)])
  val <- -get("lambda0",penden.env)*my.positive.definite.solve(get("Derv2.pen",penden.env))%*%get("Dm",penden.env)%*%beta.val.help
  return(val)
}