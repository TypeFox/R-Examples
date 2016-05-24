"find.beta" <-
function(m,alpha)
  {  
    quant.ecd <- -log(0.5*log(1/(1-alpha))) 
    c.n <- 2*log(log(m))+0.5*log(log(log(m)))-0.5*log(4*pi)
    b.n <- sqrt(2*log(log(m)))
    beta <- 1/sqrt(m)*(quant.ecd+c.n)/b.n 
  }

