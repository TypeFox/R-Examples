mantelPower <- 
function(obj, effect.size = seq(0, 1, length.out = 50), 
	alpha = 0.05)
{
    if(!inherits(obj, "mantelTest"))
	 stop("'obj' must be an object of class 'mantelTest'")
    if(alpha < 0 || alpha > 1)
	 stop("'alpha' must a numeric object between 0 and 1!")

    # Power
    nullcor <- obj$nullcor
    n <- length(nullcor)
    power <- c()
    if (obj$alternative == "greater") {
       qu <- quantile(nullcor, p = (1 - alpha)) # quantile
       for(i in 1:length(effect.size)) {
	    power[i] <- mean(nullcor >= qu - effect.size[i])
       }
    } else if (obj$alternative == "two.sided") {
       qu1 <- quantile(nullcor, p = alpha) # quantile1   
       qu2 <- quantile(nullcor, p = 1 - alpha) # quantile2 
       for(i in 1:length(effect.size)) {
          if (effect.size[i] >= 0) {
             power[i] <- mean(nullcor >= qu2 - effect.size[i])
          } else {
             power[i] <- mean(nullcor <= qu1 - effect.size[i])
          }
       }
    } else {
       qu <- quantile(nullcor, p = alpha) # quantile
       for(i in 1:length(effect.size)) {
	    power[i] <- mean(nullcor <= qu - effect.size[i])
       }
    }
    
    # Output
    out <- data.frame(effect.size = effect.size, 
       power = power)
    return(out)
}
