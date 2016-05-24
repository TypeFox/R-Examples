asypow.noncent <- function(theta.ha, info.mat, constraints, nobs.ell=1,
                              get.ho=TRUE) {
#----------------------------------------------------------------------
#
# theta.ha: Array of parameter values under the alternative hypothesis.
#
# info.mat: The information matrix, the second derivate matrix of the
#           expected log likelihood under the alternative hypothesis.
#           The negative of the hessian matrix.
#
# constraints: the constraints which set the null hypothesis from the
#     alternative hypothesis. They are in matrix form.
#          CONSTRAINT[,1] is 1 for setting parameter to a value
#                            2 for equality of two parameters
#          CONSTRAINT[,2] is case on CONSTRAINT[,1]
#               (1) Index of parameter to set to value
#               (2) Index of one of two parameters to be set equal
#          CONSTRAINT[,3] is case on CONSTRAINT[,1]
#               (1) Value to which parameter is set
#               (2) Index of other of two parameters to be set equal
#
# nobs.ell: The number of observations used in computing the information
#           matrix.  That is, info.mat is that for nobs.ell observations.
#
# get.ho: If TRUE, estimates of the parameter values under the null
#      hypothesis are calculated and returned, otherwise not.
#
#
# RETURNS a list of:
#
# w: noncentrality parameter for 1 observation
#
# df: degrees of freedom of the test
#
# theta.ho: estimates of the parameter values under the null hypothesis
#
#----------------------------------------------------------------------

      if (is.matrix(theta.ha)) theta.ha <- as.vector(t(theta.ha))
      p <- length(theta.ha)

      if (length(info.mat)==1) info.mat <- matrix(info.mat, nrow=1, ncol=1)
      dimi <- dim(info.mat)
      if (dimi[1] != dimi[2]) stop("info.mat is not a square matrix")
      if (p != dimi[1])
           stop("length of theta.ha must match dimension of info.mat")

      Aans <- asypow.construct.a(constraints,length(theta.ha))
      Ainv <- solve(Aans$a)
      hess <- t(Ainv) %*% info.mat %*% Ainv

      constraints <- Aans$constraints
      s <- dim(constraints)[1]
      r <- p - s
      
      hess.phiphi <- hess[1:s, 1:s]
      if (r > 0) {
           hess.philam <- hess[1:s, (s+1):p]
           hess.lamphi <- hess[(s+1):p, 1:s]
           hess.lamlam <- hess[(s+1):p, (s+1):p]

           w <- hess.phiphi - 
                          hess.philam %*% solve(hess.lamlam) %*% hess.lamphi
      }      else w <- hess.phiphi

      phi.ha.tran <- transformPhi(theta.ha, constraints)

      phi.diff <- Aans$phi.ho - phi.ha.tran

      w <- t(phi.diff) %*% w %*% phi.diff

#  now convert to the w statistic for one observation
      w <- w / nobs.ell

      df <- s

      if(r > 0 & get.ho) {
           lam.ho.tran <- theta.ha[-Aans$ix.con] + 
                           solve(hess.lamlam) %*% hess.lamphi %*% phi.diff

           theta.ho <- asypow.theta.ho(lam.ho.tran, constraints)

           return(list(w=w, df=df, theta.ho=theta.ho))
      }      else return(list(w=w, df=df))
}
