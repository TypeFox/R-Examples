# Added 3/28/14 for asymptotic method
medciAsymp <-
  function(mu.x,mu.y,se.x,se.y,rho=0,alpha=.05,...){
    p <- alpha/2
    quantMean <- mu.x*mu.y + se.x*se.y*rho 
    quantSE <- sqrt(se.y^2*mu.x^2+se.x^2*mu.y^2+2*mu.x*mu.y*rho*se.x*se.y+se.x^2*se.y^2+se.x^2*se.y^2*rho^2); 
    CI <- quantMean + c(qnorm(alpha/2),qnorm(1-alpha/2))*quantSE;
    res <- list(CI,quantMean,quantSE) 
    names(res) <- c(paste((1-alpha/2)*100,"% ","CI",sep=""), "Estimate", "SE")
    return(res)
  }