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

neldermead.scaletocenter <- function(this=NULL,simplex0=NULL,x0=NULL){

  hasnlcons <- optimbase.hasnlcons(this=this$optbase)
  nbve <- optimsimplex.getnbve(this=simplex0)
  xref <- transpose(optimsimplex.getx(this=simplex0,ive=1))
  for (ive in 2:nbve){
    xref <- transpose(optimsimplex.xbar(this=simplex0,iexcl=ive:nbve))
    # Transpose, because optimsimplex returns row vectors
    x <- transpose(optimsimplex.getx(this=simplex0,ive=ive))
    this <- neldermead.log(this=this,
                           msg=sprintf(paste('Scaling vertex #%d/%d at [',strvec(x),']... ',sep=''),
                                       ive,nbve))
    # Transpose x into a row vector
    tmp <- scaleinconstraints(this=this,x=x,xref=xref)
      this <- tmp$this
      status <- tmp$isscaled
      xp <- tmp$p
    rm(tmp)

    if (!status)
      stop(sprintf('neldermead.startup: Impossible to scale the vertex #%d/%d at [%s] into inequality constraints',
                   ive,nbve,strvec(x)),
           call.=FALSE)
    if (any(x!=xp)){
      tmp <- optimbase.function(this=this$optbase,x=xp,index=2)
      if (hasnlcons){
        fv <- tmp$f
        c <- tmp$c
        index <- tmp$index
      } else {
        fv <- tmp$f
        index <- tmp$index
      }
      # Transpose xp, which is a column vector
      simplex0 <- optimsimplex.setve(this=simplex0,ive=ive,fv=fv,x=transpose(xp))
    }
  }
  return(simplex0)
}

