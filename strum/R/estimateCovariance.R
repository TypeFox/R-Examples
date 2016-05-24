#==============================================================================
# File: estimateCovariance.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver to calculate the covariance matrix of parameters.
#        -  Calculate covariance matrix of Delta parameters
#        -  Calculate covariance matrix of Theta parameters
#        -  Calculate covariance matrix of differance = (delta - theta)
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Estimate Covariance(delta) = (J'FJ'), Equation 8.
#------------------------------------------------------------------------------
.estimateDeltaCovariance = function(deltaDeriv)
{
  #
  # 1. Estimate J and F
  #---------------------
  K = ncol(deltaDeriv$fDeriv) # equal to numPed

  # Use .bdiag to stack them on a diagonal
  # - a sparse matrix is constructed.
  # Needs library(Matrix)
  J = .bdiag(deltaDeriv$sDeriv)

  if( deltaDeriv$numTrt > 1 )
  {
    index = 1
    rowOffset = deltaDeriv$numCov + deltaDeriv$numVC

    for( t1 in 2:deltaDeriv$numTrt )
    {
      r1e = rowOffset * t1 + .getCovBlockCount(t1-2) * deltaDeriv$numVC
      r1s = r1e - rowOffset + 1

      for( t2 in (t1-1):1 )
      {
        r2e = rowOffset * t2 + .getCovBlockCount(t2-2) * deltaDeriv$numVC
        r2s = r2e - rowOffset + 1

        cs = r1e + 1 + (t1-t2-1) * deltaDeriv$numVC 
        ce = cs + deltaDeriv$numVC - 1

        # upper matrix
        #J[r1s:r1e, cs:ce] = deltaDeriv$s1Deriv[[index]]
        #J[r2s:r2e, cs:ce] = deltaDeriv$s2Deriv[[index]]

        # lower matrix
        J[cs:ce, r1s:r1e] = t(deltaDeriv$s1Deriv[[index]])
        J[cs:ce, r2s:r2e] = t(deltaDeriv$s2Deriv[[index]])

        index = index + 1
      }
    }
  }

  J = J/K

  F = (deltaDeriv$fDeriv %*% t(deltaDeriv$fDeriv))/K

  #
  # 2. Calculate covariance matrix
  #--------------------------------
  #Ji = solve(J)
  Ji = tryCatch(solve(J),
                error=function(e)
                      {
                        stop(paste("Covariance matrix of the model parameters ",
                                   "is not invertable in step 1!  The model may ",
                                   "not be identifiable, or the sample size ",
                                   "is too small.",
                             sep=""))
                      })

  JiFJit = Ji %*% F %*% t(Ji)

  w = 1/diag(JiFJit)

  return(list(K=K, deltaCov=JiFJit, w=w))
}

#------------------------------------------------------------------------------
# Estimate asymptotic covariance(theta) = (A'BA')/k of the model.
#------------------------------------------------------------------------------
.estimateThetaCovariance = function(model, s1, theta)
{
  dFdTheta = .funDeriv(fun=function(theta) .thetaToDelta(model,theta), theta)
  #dfOrthog = Null(t(dFdTheta)) #Library MASS
  dfOrthog = .strumNull(t(dFdTheta), tol=1e-07)

  # Calculate covariance matrix of theta.
  #---------------------------------------
  W = diag(s1$w)
  WJFJW = W %*% s1$deltaCov %*% W 

  A = dFdTheta %*% W     %*% t(dFdTheta)
  B = dFdTheta %*% WJFJW %*% t(dFdTheta)

  #Ai = solve(A)
  Ai = tryCatch(solve(A), error=function(e) return(ginv(A)))

  thetaCov = (Ai %*% B %*% Ai) / s1$K

  # Calculate covariance matrix of diff=(delta-theta).
  #----------------------------------------------------
  diff     = s1$delta - .thetaToDelta(model, theta)
  dFAidFW  = t(dFdTheta) %*% Ai %*% dFdTheta %*% W
  IdFAidFW = diag(nrow(W)) - dFAidFW

  diffCov = (IdFAidFW %*% s1$deltaCov %*% t(IdFAidFW))

  U = W - (W %*% dFAidFW)

  return(list(thetaCov = thetaCov,
              diff     = diff,
              diffCov  = diffCov,
              dfOrthog = dfOrthog,
              U        = U))
}

#------------------------------------------------------------------------------
# Find the number of previous covariance blocks.
#------------------------------------------------------------------------------
.getCovBlockCount = function(N)
{
  if( N <= 0 )
    return(0)
  else
    return(sum(1:N))
}
