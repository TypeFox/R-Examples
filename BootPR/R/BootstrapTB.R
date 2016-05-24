BootstrapTB <-
function(x,p,h,nboot)
{
    #set.seed(12345)
    B <- LSMBT(x,p)
    n <- nrow(x)
    b <- B$coef 
    e <- sqrt( (n-p) / ( (n-p)-length(b)))*B$resid
    btem1 <- numeric(length(b))    
    for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(e)))
        es <- e[index,1]
        xs <- ysbT(x, b, es)
        btem1 <- btem1 + LSMBT(xs,p)$coef/nboot
    }
    bc <- 2*b-btem1
    
    bc <- adjust(b,bc,p)
    
    if(sum(b) != sum(bc))
    bc[(p+1):(p+2),] <- RE.LSMTB(x,p,bc)
    
    e <- RESIDT(x,bc)

    if(h > 0)
    f <- ART.Fore(x,bc,h)
return(list(coef=bc,resid=e,forecast=f))
}
