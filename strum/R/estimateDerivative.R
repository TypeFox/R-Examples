#==============================================================================
# File: computeDerivative.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver to calculate the derivatives of C and V at the optimum
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Estimate the 1st and 2nd derivatives of C and V
#------------------------------------------------------------------------------
.estimateDeltaDerivative = function(y, x, vc, betaMat, sigmaMat)
{
  numTrt = ncol(y$yAll[[1]])
  numCov = ncol(x$xAll[[1]])
  numVC  = length(vc$vcAll[[1]])
  numPed = length(y$yAll)
  
  fDeriv = matrix(ncol = numPed, nrow = 0)
  sDeriv = list()
  s1Deriv = list()
  s2Deriv = list()

  for( t1 in 1:numTrt )
  {	
    # 1. univariate derivatives
    #---------------------------

    bet1 = matrix(betaMat[,t1], ncol = 1)
    sig1 = lapply(sigmaMat, function(s) return(matrix(s[t1,t1])))

    # all
    #-----
    t1dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t1)
    yt1  = t1dat$y
    xt1  = t1dat$x
    vct1 = t1dat$vc

    dUni = .Call("computeDeriv",
                 yt1, xt1, vct1, list(), list(), list(), list(), bet1, sig1)

    if( length(y$yPro) )
    {
      # probands
      #----------
      pt1dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t1)
      ypt1  = pt1dat$y
      xpt1  = pt1dat$x
      vcpt1 = pt1dat$vc

      dUniP = .Call("computeDeriv",
                    ypt1, xpt1, vcpt1, list(), list(), list(), list(), bet1, sig1)

      fDeriv = rbind(fDeriv, dUni$dPar - dUniP$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par - dUniP$d2Par
    } else
    {
      fDeriv = rbind(fDeriv, dUni$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par
    }

    if( t1 == 1 )
      next

    # 2. bivariate derivatives
    #--------------------------
    for( t2 in (t1-1):1 )
    {
      bet12 = matrix(betaMat[,c(t1,t2)], ncol=2)
      sig12 = lapply(sigmaMat, function(s) return(s[c(t1,t2),c(t1,t2)]))

      # all
      #-----
      t2dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t2)
      yt2  = t2dat$y
      xt2  = t2dat$x
      vct2 = t2dat$vc
      vc12 = .getVCFilteredData(vc$vcAll, t1dat$missing, t2dat$missing)

      dBi = .Call("computeDeriv",
                  yt1, xt1, vct1, yt2, xt2, vct2, vc12, bet12, sig12)

      if( length(y$yPro) )
      {
        # probands
        #----------
        pt2dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t2)
        ypt2  = pt2dat$y
        xpt2  = pt2dat$x
        vcpt2 = pt2dat$vc
        vcp12 = .getVCFilteredData(vc$vcPro, pt1dat$missing, pt2dat$missing)
 
        dBiP = .Call("computeDeriv",
                     ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12, bet12, sig12)

        fDeriv = rbind(fDeriv, dBi$dPar - dBiP$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par - dBiP$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS - dBiP$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS - dBiP$dP2dS
      } else
      {
        fDeriv = rbind(fDeriv, dBi$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS
      }      
    }
  }

  return(list(numTrt=numTrt, numCov=numCov, numVC=numVC,
              fDeriv=fDeriv, sDeriv=sDeriv, s1Deriv=s1Deriv, s2Deriv=s2Deriv))
}
