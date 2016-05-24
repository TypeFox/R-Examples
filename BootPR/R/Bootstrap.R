Bootstrap <-
function(x,p,h,nboot)
{
    #set.seed(12345)
    B <- LSM(x,p)
    n <- nrow(x)
    b <- B$coef 
    e <- sqrt( (n-p) / ( (n-p)-length(b)))*B$resid
    btem1 <- numeric(length(b))    
    for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(e)))
        es <- e[index,1]
        xs <- ys(x, b, es)
        btem1 <- btem1 + OLS.AR(xs,p,0,0)$coef/nboot
    }
    bc <- 2*b-btem1
   
    bc <- adjust(b,bc,p)
      
    if(sum(b) != sum(bc))
    bc[p+1] <- mean(x)*(1-sum(bc[1:p]))
    
    e <- RESID(x,bc)
    f <- {}
    if(h > 0)
    f <- AR.Fore(x,bc,h)
return(list(coef=bc,resid=e,forecast=f))
}
