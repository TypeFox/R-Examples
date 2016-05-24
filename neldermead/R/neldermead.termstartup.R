# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the neldermead component
# originally written by Michael Baudin for Scilab :
# "Nelder-Mead User's Manual", 2010, Consortium Scilab - Digiteo,
# Michael Baudin, http://wiki.scilab.org/The_Nelder-Mead_Component

neldermead.termstartup <- function(this=NULL){
  #
  # Criteria #8 : Kelley stagnation, based on
  # a sufficient decrease condition
  #
  if (this$kelleystagnationflag){
    if (!this$kelleynormalizationflag){
      this$kelleyalpha <- this$kelleystagnationalpha0
    } else {
      sg <- optimsimplex.gradientfv(this=this$simplex0,fun=neldermead.costf,method='forward',data=this)
      nsg <- transpose(sg) %*% sg
      sigma0 <- optimsimplex.size(this=this$simplex0,method='sigmaplus')
      if (nsg==0.0){
        this$kelleyalpha <- this$kelleystagnationalpha0
      } else {
        this$kelleyalpha <- this$kelleystagnationalpha0 * sigma0 / nsg
      }
    }
    this <- neldermead.log(this=this,
                           msg=sprintf('Test Stagnation Kelley : alpha0 = %e',this$kelleyalpha))
  }
  #
  # Criteria #10 : variance of function values
  #
  if (this$tolvarianceflag)
    this$variancesimplex0 <- optimsimplex.fvvariance(this=this$simplex0)

  return(this)

}

