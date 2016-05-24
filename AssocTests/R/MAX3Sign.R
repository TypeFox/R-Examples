#Calculate statistics
max3Sign <- function(T.max, cor.rec.add, cor.rec.dom, cor.add.dom)
{          
    sigma <- matrix(c(1,cor.rec.dom,cor.rec.dom,1),ncol=2) 
    
    let <- 1-cor.rec.dom^2
    w0 <- (cor.rec.add - cor.rec.dom * cor.add.dom)/let
    w1 <- (cor.add.dom - cor.rec.dom * cor.rec.add)/let
    
    wet <- sqrt(let)
    f1 <- function(z0)
    {
        L <- length(z0)
        qet <- double(L)
        for (i in 1:L)
        {
            qet[i] <- pnorm((T.max-cor.rec.dom*z0[i])/wet) * dnorm(z0[i])
        }
        
        qet
    }
    
    f2 <- function(z0)
    {
        L <- length(z0)
        qet <- double(L)
        for (i in 1:L)
        {
            qet[i] <- pnorm(((T.max-w0*z0[i])/w1-cor.rec.dom*z0[i])/wet) * dnorm(z0[i])
        }
        
        qet
    }
    
    f3 <- function(z0)
    {
        L <- length(z0)
        qet <- double(L)
        for (i in 1:L)
        {
            qet[i] <- pnorm((-T.max-cor.rec.dom*z0[i])/wet) * dnorm(z0[i])
        }
        
        qet
    }
    
    x1 <- integrate(f1, 0, T.max*(1-w1)/w0, subdivisions=10000, rel.tol=1e-10)$value
    x2 <- integrate(f2, T.max*(1-w1)/w0, T.max,subdivisions=10000, rel.tol=1e-10)$value
    x3 <- integrate(f3, 0, T.max, subdivisions=10000, rel.tol=1e-10)$value
    
    p.value <- 1-2*(x1+x2-x3)
    
    p.value    
}
