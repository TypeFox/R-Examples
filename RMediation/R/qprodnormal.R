qprodnormal <-
function(p, mu.x, mu.y, se.x, se.y, rho=0, lower.tail=TRUE, type="dop", n.mc=1e5){
  if(p<=-1 || p>=1)
    stop("p must be between -1 and 1!")
  if(!is.numeric(mu.x))
    stop("Argument mu.x must be numeric!")
  if(!is.numeric(mu.y))
    stop("Argument mu.y must be numeric!")
  if(!is.numeric(se.x))
    stop("Argument se.x must be numeric!")
  if(!is.numeric(se.y))
      stop("Argument se.y must be numeric!")
  if(!is.numeric(rho))
      stop("Argument rho  must be numeric!")
  if(rho<=-1 || rho>=1)
    stop("rho must be between -1 and 1!")
  if(!is.numeric(n.mc) || is.null(n.mc))
    n.mc=1e5 # sets n.mc to default

  if (type=="all" || type=="All" || type=="ALL")
    {
       ##cat("Meeker method:\n")
        q2 <- qprodnormalMeeker(p, mu.x, mu.y, se.x, se.y, rho, lower.tail)
        ##cat("Monte Carlo method:\n")
        q3 <- qprodnormalMC(p, mu.x, mu.y, se.x, se.y, rho, lower.tail, n.mc)
        res <- list(q2,q3)
        names(res) <- c( "Distribution of Product", "Monte Carlo")
        return(res)
    }
  else if (type=="DOP" || type=="dop")
    {
        ##cat("Meeker method:\n")
        q2 <- qprodnormalMeeker(p, mu.x, mu.y, se.x, se.y, rho, lower.tail)
        return(q2)
    }
  else if (type=="MC" || type=="mc" || type=="Mc")
    {
        ##cat("Monte Carlo method:\n")
        q3 <- qprodnormalMC(p, mu.x, mu.y, se.x, se.y, rho, lower.tail, n.mc)
        return(q3)
    }
  else stop("Wrong type! please specify type=\"all\", \"DOP\", or \"MC\" ")
}


