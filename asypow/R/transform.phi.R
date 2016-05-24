transformPhi <- function(theta.ha, constraints) {

      s <- dim(constraints)[1]
      phi.ha.tran <- rep(NA,s)

      for(i in 1:s) {
           if(constraints[i,1] == 1)
                phi.ha.tran[i] <- theta.ha[constraints[i,2]] else 
                   phi.ha.tran[i] <- 
                     theta.ha[constraints[i,3]] - theta.ha[constraints[i,2]]
      }

      return(phi.ha.tran)
}
