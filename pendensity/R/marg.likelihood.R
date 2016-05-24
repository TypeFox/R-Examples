marg.likelihood <- function(penden.env,pen.likelihood) {
  help <- eigen(get("Dm",penden.env),symmetric=TRUE)

  index <- which(help$values>1e-5)
  evalues <- help$values[index]
  Utilde <- help$vectors[,index]

  k1 <- 0.5*sum(log(get("lambda0",penden.env)*evalues))
  k2 <- pen.likelihood

  eneu <- eigen(t(Utilde)%*%get("Derv2.pen",penden.env)%*%Utilde)$values

  if(any(eneu<=0)) print("eneu <=0")
  k3 <- -0.5*sum(log(eneu))
  return(k1+k2+k3)
}
