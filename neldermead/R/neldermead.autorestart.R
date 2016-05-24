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

neldermead.autorestart <- function(this=NULL){

  restartmax <- this$restartmax
  reached <- FALSE
  for (iloop in 1:(restartmax+1)){
    this <- neldermead.log(this=this,
                           msg='*****************************************************************')
    this <- neldermead.log(this=this,msg=sprintf('Try #%d/%d.',iloop,restartmax+1))
    #
    # Run algorithm
    this <- neldermead.algo(this=this)

    #
    # Must we restart ?
    tmp <- neldermead.istorestart(this=this)
      this <- tmp$this
      istorestart <- tmp$istorestart
    rm(tmp)
    if (istorestart){
      this <- neldermead.log(this=this,'Must restart.')
    } else {
      this <- neldermead.log(this=this,'Must not restart.')
    }
    if (!istorestart){
      reached <- TRUE
      break
    }
    if (iloop<restartmax+1){
      # We are going to perform a restart
      this$restartnb <- this$restartnb + 1
      this <- neldermead.log(this=this,msg='Updating simplex.')
      this <- neldermead.updatesimp(this=this)
    }
  }
  if (reached){
    this <- neldermead.log(this=this,
                           msg=sprintf('Convergence reached after %d restarts.',this$restartnb))
  } else {
    this <- neldermead.log(this=this,
                           msg=sprintf('Convergence not reached after maximum %d restarts.',this$restartnb))
    this$optbase <- optimbase.set(this=this$optbase,key='status',value='maxrestart')
  }
  return(this)
}

