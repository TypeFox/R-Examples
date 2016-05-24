"composite2COP" <-
function(u,v,para,...) {
  alpha <- para$alpha; alpha.p <- 1 - alpha
  beta  <- para$beta;   beta.p <- 1 - beta
  return(COP(cop=para$cop1, u^alpha,   v^beta,   para=para$para1) *
         COP(cop=para$cop2, u^alpha.p, v^beta.p, para=para$para2))
}

