## ' Computes the Expected reduction of the volume of the excursion sets behind a 2D Pareto front.
## ' @title Expected reduction of the volume of excursion
## ' @param phi.x.bar see below,
## ' @param phi.x.tilde see below,
## ' @param  phi.eta.x see below, 
## ' @param  phi.y.bar see below, 
## ' @param  phi.y.tilde see below, 
## ' @param  phi.eta.y see below, 
## ' @param  phi2.x.x see below, 
## ' @param  phi2.x.eta see below, 
## ' @param  phi2.y.y see below, 
## ' @param  phi2.y.eta vectors of size n.pareto*n.integration.points. 'x' refers to the first objective, 'y' to the second.
## ' @details To be called by \code{\link[GPareto]{crit_SUR}}
## ' 
## ' @return list with three elements: 
## ' \code{crit} is the expected reduction
## ' \code{pijold} is the current probability of non-domination
## ' \code{pij} is the expected future probability of non-domination
## ' @references
## ' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
## ' \emph{Statistics and Computing}
## ' @export

EEV.2D.computation <- function(phi.x.bar, phi.x.tilde, phi.eta.x, phi.y.bar, phi.y.tilde, phi.eta.y, 
                               phi2.x.x, phi2.x.eta, phi2.y.y, phi2.y.eta)
{
  ############################################################################
  n.pareto <- length(phi.x.bar)
  n.integration.points <- ncol(phi.x.tilde)
  
  #*** Current excursion volume **********************************************
  oldcrit <- crit <- 0
  pij <- pijold  <- matrix(0, n.pareto+1, n.integration.points)
  
  pijold[1,] <- phi.x.tilde[1,]
  if (n.pareto > 1){
    for (j in 2:(n.pareto)){
      pijold[j,] <- (phi.x.tilde[j,] - phi.x.tilde[j-1,])*phi.y.tilde[j-1,]
    }
  }
  pijold[n.pareto+1,] <- (1 - phi.x.tilde[n.pareto,])*phi.y.tilde[n.pareto,]
  
  #*** New excursion volume **************************************************
  # j = 1
  phipp.obj1   <- phi2.x.x[1,]
  phipnu.obj1  <- phi2.x.eta[1,]
  phi.ftilde.p <- phi.x.tilde[1,]
  
  pij[1,] <- ((phipp.obj1 - phipnu.obj1)*(phi.eta.y - 1) + phi.ftilde.p)
  
  # j = 2 ... n.pareto
  if (n.pareto > 1){
    for (j in 2:(n.pareto)){
      phimm.obj1  <- phipp.obj1
      phimnu.obj1 <- phipnu.obj1
      phi.ftilde <- phi.ftilde.p
      
      phipp.obj1  <- phi2.x.x[j,]     
      phipnu.obj1 <- phi2.x.eta[j,]
      phipp.obj2  <- phi2.y.y[j-1,] 
      phipnu.obj2 <- phi2.y.eta[j-1,]
      
      phi.gtilde   <- phi.y.tilde[j-1,]
      phi.ftilde.p <- phi.x.tilde[j,]
      
      pij[j,] <-  ((phipp.obj1 - phipnu.obj1 + phimnu.obj1 - phimm.obj1)*(phipnu.obj2 - phipp.obj2) + (phi.ftilde.p - phi.ftilde)*phi.gtilde )
    }
  }
  # j = n.pareto+1
  j = n.pareto
  phimnu.obj1  <- phipnu.obj1
  phimm.obj1   <- phipp.obj1
  phipp.obj2  <- phi2.y.y[j,]
  phipnu.obj2 <-  phi2.y.eta[j,]
  phi.ftilde <- phi.ftilde.p
  phi.gtilde <- phi.y.tilde[j,]
  
  pij[j+1,] <- ((1 - phi.eta.x + phimnu.obj1 - phimm.obj1)*(phipnu.obj2 - phipp.obj2) + (1 - phi.ftilde)*phi.gtilde )
  
  crit <- sum(pijold - pij)
  return(list(crit, pijold, pij))
}