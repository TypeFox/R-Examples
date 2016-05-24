"convex2COP" <-
function(u,v, para, ...) {
  return(   para$alpha  * COP(u, v, cop=para$cop1, para=para$para1) +
         (1-para$alpha) * COP(u, v, cop=para$cop2, para=para$para2))
}

