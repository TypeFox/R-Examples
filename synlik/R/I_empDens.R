# Calculating log-likelihood
.empDens <- function(y, X, saddle, decay, tol = 1e-6, log = TRUE)
{
  
#   # If saddle == TRUE returns saddlepoint density, otherwise a normal density
#   if(saddle == TRUE)
#   {
#     
#     tmp <- dsaddle(y = y, X = X, tol = tol, decay = decay, log = TRUE)
#     llk <- tmp$llk
#     mix <- tmp$mix
#     
#   } else  {
    
    tmp <- robCov( t(X) )
    
    # If there are some statistics with zero variace we remove them
    if( length(tmp$lowVar) ) y <- y[-tmp$lowVar]
    
    rss <- sum( (tmp$E%*%as.vector(y-tmp$mY))^2 )
    llk <- -rss/2 - tmp$half.ldet.V - 0.5 * length(y) * log(2 * pi)
    mix <- 0
    
#   }
  
  list("logLik" = ifelse(log, llk, exp(llk)), "mix" = mix)
}
