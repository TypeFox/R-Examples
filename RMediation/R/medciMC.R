medciMC <-
  function(mu.x, mu.y, se.x, se.y, rho=0, alpha=.05,  n.mc=1e7){
    mean.v <- c(mu.x,mu.y)
    var.mat <- matrix(c(se.x^2,se.x*se.y*rho,se.x*se.y*rho,se.y^2),2)
    a_b <- matrix(rnorm(2*n.mc),ncol=n.mc)
    a_b <- crossprod(chol(var.mat),a_b)+mean.v
    a_b <- t(a_b)
    ab <- a_b[,1]*a_b[,2]
    quantMean <- mean(ab)
    quantSE <- sd(ab)
    quantError <- quantSE/n.mc
    CI <- (quantile(ab,c(alpha/2,1-alpha/2))) # Added 3/28/14
    #names(CI) <- c(paste((alpha/2*100),"%"),paste((1-alpha/2)*100,"%"))
    res <- list(CI,quantMean,quantSE,quantError) # Added 3/28/14
    names(res) <- c(paste((1-alpha/2)*100,"% ","CI",sep=""), "Estimate", "SE","MC Error")
    return(res)
  }
