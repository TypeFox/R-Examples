info.poisson.design <- function(model="linear", theta, xpoints,
                             natx=1, group.size=1) {
#-----------------------------------------------------------------------
#        Returns the information matrix for a poisson design
#
# model: One of {"linear", "quadratic"} Only enough to ensure a unique
#           match need be supplied.
#
# theta: Matrix of parameters of the linear part of the model.
#             Each row represents a group.
#
# xpoints: Matrix of covariate values for each group.
#             Each row represents a group.
#
# natx: Matrix of number of observations at xpoints for each group.
#           Each row represents a group.
#           At covariate value xpoint[i,j] there are natx[i,j] observations.
#
# group.size: The relative number of observations in each group
#
#
# Returns: The the information matrix for one observation for this design.
#     The observation is assumed to be spread over xpoints in proportion
#     to natx.
#
#-----------------------------------------------------------------------

      if (is.vector(theta)) theta <- t(as.matrix(theta))
      dimt    <- dim(theta)
      ngroups <- dimt[1]

      if (is.vector(xpoints)) {
           xpoints <- matrix(xpoints,ngroups,length(xpoints), byrow=TRUE)
           dimp <- dim(xpoints)
      } else {
           dimp <- dim(xpoints)
           if (dimp[1] != ngroups)
                stop("Number of rows of xpoints and theta must match")
      }
      if (is.vector(natx)) {
           if (length(natx) == 1)
                natx <- matrix(natx,ngroups,dimp[2],byrow=TRUE) else {
                if (length(natx) != dimp[2] )
                     stop ("length of natx must match number of xpoints")
                natx <- matrix(natx,ngroups,dimp[2],byrow=TRUE)
           }
           dimn <- dim(natx)
      } else {
           dimn <- dim(natx)
           if (dimn[1] != dimp[1] || dimn[2] != dimp[2])
                stop("xpoints and natx must have the same dimensions")
      }

      lngrpsz <- length(group.size)
      if (lngrpsz == 1) group.size <- rep(group.size,ngroups) else 
            if (ngroups != lngrpsz)
       stop("\nNumber of rows of theta and length of group.size must match")

      info <- vector("list",dimp[1])

      for (j in 1:dimp[1]) {
           info[[j]] <- 0
           for (i in 1:length(xpoints[j,]))
                if(natx[j, i] != 0) info[[j]] <- info[[j]] + natx[j,i] *
                      info.poisson.one(model,theta[j,],xpoints[j,i])
           info[[j]] <- info[[j]] * group.size[j]
      }

      info <- k.blocks.info(info)

      return(info/sum(group.size * natx))
}
