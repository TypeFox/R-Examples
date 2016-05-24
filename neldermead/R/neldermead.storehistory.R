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

neldermead.storehistory <- function(this=NULL,n=NULL,fopt=NULL,xopt=NULL,
                                    fv=NULL,xcoords=NULL){
  
  storehistory <- optimbase.get(this=this$optbase,key='storehistory')
  iterations <- optimbase.get(this=this$optbase,key='iterations')
  verbose <- neldermead.get(this=this,key='verbose')
  nbve <- neldermead.get(this=this,key='simplex0')$nbve
  if (storehistory){
    this$optbase <- optimbase.histset(this=this$optbase,
                                      iter=iterations,
                                      key='historyfopt',
                                      value=fopt)
    this$optbase <- optimbase.histset(this=this$optbase,
                                      iter=iterations,
                                      key='historyxopt',
                                      value=xopt[1:n])
    this$historysimplex[[iterations]] <- simplex(verbose=verbose,
      x=xcoords[1:nbve,1:n,drop=FALSE],
      fv=fv,
      n=n,
      nbve=nbve)
  }
  return(this)
}

