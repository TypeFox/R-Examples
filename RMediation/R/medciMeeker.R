medciMeeker <-
function(mu.x,mu.y,se.x,se.y,rho=0,alpha=.05,...){
  p <- alpha/2
  q.l <- qprodnormalMeeker(p, mu.x=mu.x, mu.y=mu.y, se.x=se.x, se.y=se.y, rho=rho, lower.tail=TRUE)
  q.u <- qprodnormalMeeker(p, mu.x=mu.x, mu.y=mu.y, se.x=se.x, se.y=se.y, rho=rho, lower.tail=FALSE)
  CI <- c(q.l[[1]],q.u[[1]]) #confidence interval
  quantMean <- mu.x*mu.y + se.x*se.y*rho 
  quantSE <- sqrt(se.y^2*mu.x^2+se.x^2*mu.y^2+2*mu.x*mu.y*rho*se.x*se.y+se.x^2*se.y^2+se.x^2*se.y^2*rho^2); 
  res <- list(CI,quantMean,quantSE) 
  names(res) <- c(paste((1-alpha/2)*100,"% ","CI",sep=""), "Estimate", "SE")
  return(res)
}

