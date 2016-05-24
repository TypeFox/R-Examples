## modify.normalize marker
ModifyNormalization <- function (x)
{
    x <-as.vector(x)
    dex <- which(!is.na(x))
    
    row.sum <- sum(x[dex])
    row.valid <- length(dex)    

    row.mean <- row.sum/row.valid
    p <- (row.sum+0.5)/(1+row.valid) 
  
    y <- rep(0,length(x))
    y[dex] <- (x[dex]-row.mean)/sqrt(p*(1-p))
            
    y
}
