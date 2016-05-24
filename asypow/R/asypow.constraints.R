asypow.constraints <- function(constraints) {
#-----------------------------------------------------------------------
#
#       Fixes the constraints, makes sure they are consistant
#
#-----------------------------------------------------------------------

      if (is.vector(constraints)) {
           if (length(constraints) != 3)
                     stop("constraint vector must be of length 3")
           constraints <- t(as.matrix(constraints))
           s <- 1
      }      else {
           dimc <- dim(constraints)
           if (dimc[2] != 3)
                     stop("constraint matrix must have 3 columns")
           s <- dimc[1]
      }

      if (any(constraints[,1] != 1 & constraints[,1] != 2))
                     stop("bad constraint matrix")

      for (i in 1:s) {
           if (constraints[i,1] == 2) {
                if (constraints[i,2] == constraints[i,3])
                          stop("bad constraint matrix") else 
                     if (constraints[i,3] < constraints[i,2])
                     constraints[i,] <- c(constraints[i,1],constraints[i,3],
                                                             constraints[i,2])
           }
      }

      if (any(duplicated(constraints[,2]))) stop("bad constraint matrix")

      ord             <- order(constraints[,2])
      constraints[,1] <- constraints[ord,1]
      constraints[,2] <- constraints[ord,2]
      constraints[,3] <- constraints[ord,3]

      min.no.cons <- max(c(constraints[,2],constraints[constraints[,1]==2,3]))

      return(list(constraints=constraints,min.no.cons=min.no.cons))

}
