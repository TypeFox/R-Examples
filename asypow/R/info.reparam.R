info.reparam <- function(theta, info.mat, dg) {
#-----------------------------------------------------------------------
#        Returns the information matrix for a binomial design
#
# theta: Matrix of parameters of the linear part of the model.
#          Each row represents a group. This is under the original
#          parameterization.
#
# info.mat: The information matrix under the original parameterization
#
# dg: A function that computes the partial deravatives of g,the
#          transformation function. Takes theta and returns dg/dtheta[1]
#          dg/dtheta[2] ...
#
#
# Returns: The information matrix under the new parameterization
#
#-----------------------------------------------------------------------

      if (is.vector(theta)) theta <- t(as.matrix(theta))
      dimt <- dim(theta)
      ngroups <- dimt[1]

      if(!is.matrix(info.mat)) stop("info.mat must be a matrix")
      dimi <- dim(info.mat)
      if(dimi[1]/2 != ngroups)
                stop("info.mat and theta do not match up")

      Dg <- vector("list",ngroups)
      for (j in 1:ngroups) Dg[[j]] <- dg(theta[j,])
      Dg <- k.blocks.info(Dg)

      Va <- solve(info.mat)

      var.new <- Dg %*% Va %*% t(Dg)

      info.new <- solve(var.new)

      return(info.new)
}
