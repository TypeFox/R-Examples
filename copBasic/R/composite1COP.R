"composite1COP" <-
function(u,v,para,...) {
  alpha <- para$alpha; alpha.p <- 1 - alpha
  beta  <- para$beta;   beta.p <- 1 - beta
  return(u^alpha * v^beta *
         COP(cop=para$cop1, u^alpha.p, v^beta.p, para=para$para1))
}

