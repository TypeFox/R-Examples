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

neldermead.search <- function(this=NULL){

  withderivatives <- optimbase.get(this=this$optbase,key='withderivatives')
  if (withderivatives)
    stop('neldermead.search: The -withderivatives option is true but all algorithms in neldermead are derivative-free.',
         call.=FALSE)
  if (!this$startupflag){
    this <- neldermead.startup(this=this)
    this$startupflag <- TRUE
  }
  
  neldermead.outputcmd(this=this,state='init',simplex=this$simplex0,step='init')
  if (this$restartflag){
    this <- neldermead.autorestart(this)
  } else {
    this <- neldermead.algo(this=this)
  }
  neldermead.outputcmd(this=this,state='done',simplex=this$simplexopt,step='done')

  return(this)

}

