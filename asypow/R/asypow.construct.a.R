asypow.construct.a <- function(constraints,p) {
#----------------------------------------------------------------------
#               Contructs the matrix A using constraints
#
# constraints: The constraints which set the null hypothesis from the
#     alternative hypothesis. They are in matrix form.
#          CONSTRAINT[,1] is 1 for setting parameter to a value
#                            2 for equality of two parameters
#          CONSTRAINT[,2] is case on CONSTRAINT[,1]
#               (1) Number of parameter to set to value
#               (2) Number of one of two parameters to be set equal
#          CONSTRAINT[,3] is case on CONSTRAINT[,1]
#               (1) Value to which parameter is set
#               (2) Number of other of two parameters to be set equal
#
#     p : The number of parameters in the model.
#
#
# RETURNS a list:
#
#     a : matrix that will multiply parameters
#
#     phi.ho : values of constrained parameters under the null model
#
#     ix.con : ixdex corresponding to phi.ho
#----------------------------------------------------------------------

      new.cons    <- asypow.constraints(constraints)
      constraints <- new.cons$constraints
      if (p < new.cons$min.no.cons)
                     stop("p is too small for the constraint matrix")

      s <- dim(constraints)[1]
      r <- p - s

      ix.con <- constraints[,2]
      if (r > 0) {
           ix.unc <- (1:p)[-ix.con]
           A <- matrix(0,ncol=p,nrow=r)
           for(i in 1:r) A[i,ix.unc[i]] <- 1
      }      else A <- NULL

      phi.ho <- rep(NA,s)
      Acon <- matrix(0,ncol=p,nrow=s)
      for(i in 1:s) {
           if(constraints[i,1] == 1) {
                Acon[i,constraints[i,2]] <- 1
                phi.ho[i] <- constraints[i,3]
           }  else {
                Acon[i,constraints[i,2]] <- 1
                Acon[i,constraints[i,3]] <- -1
                phi.ho[i] <- 0
           }
      }

      return(list(a=rbind(Acon,A),phi.ho=phi.ho,ix.con=ix.con,
                                              constraints=constraints))
}
