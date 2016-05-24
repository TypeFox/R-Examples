"fx" <-
function(bta,n,alpha,p1,s)
 {
 pb <- p1 - s*bta
 f <- 2*n*(pb*log(pb) + (1-pb)*log(1-pb) - log(0.5)) - (qnorm(alpha))^2 
 return(f)
 }

