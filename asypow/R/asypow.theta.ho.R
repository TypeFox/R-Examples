asypow.theta.ho <- function(lam.ho, constraints) {

      p <- length(lam.ho) + length(constraints[,1])
      theta.ho <- rep(NA,p)

      ix.con <- constraints[,2]
      theta.ho[-ix.con] <- lam.ho

      if(any(constraints[,1]==1)) {
           ix.set <- constraints[constraints[,1]==1,2]
           theta.ho[ix.set] <- constraints[constraints[,1]==1,3]
      }

      if(any(constraints[,1]==2)) {
           ix.equ <- constraints[constraints[,1]==2,2]
           for(i in rev(ix.equ)) {
                theta.ho[i] <- theta.ho[constraints[i,3]]
           }
      }

      return(theta.ho)
}
