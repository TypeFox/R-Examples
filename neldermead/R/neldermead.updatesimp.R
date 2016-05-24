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

neldermead.updatesimp <- function(this=NULL){

  simplex0 <- optimsimplex()$newobj
  xopt <- optimbase.get(this=this$optbase,key='xopt')
  
  if (!any(this$restartsimplexmethod==c('axes','spendley','pfeffer','randbounds','oriented')))
    stop(sprintf('neldermead_updatesimp: Unexpected key %s',this$restartsimplexmethod),
         call.=FALSE)
  
  if (this$restartsimplexmethod=='oriented'){
    tmp <- optimsimplex(method='oriented',simplex0=this$simplexopt,fun=costf.transposex,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$restartsimplexmethod=='axes'){
    simplex0 <- optimsimplex.destroy(this=simplex0)
    tmp <- optimsimplex(method='axes',x0=transpose(xopt),fun=costf.transposex,
                            len=this$simplex0length,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$restartsimplexmethod=='spendley'){
    simplex0 <- optimsimplex.destroy(this=simplex0)
    tmp <- optimsimplex(method='spendley',x0=transpose(xopt),fun=costf.transposex,
                            len=this$simplex0length,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$restartsimplexmethod=='pfeffer'){
    simplex0 <- optimsimplex.destroy(this=simplex0)
    tmp <- optimsimplex(method='pfeffer',x0=transpose(xopt),fun=costf.transposex,
                            deltausual=this$simplex0deltausual,deltazero=this$simplex0deltazero,
                            data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$restartsimplexmethod=='randbounds'){
    if (this$boxnbpoints=='2n'){
      this$boxnbpointseff <- 2 * this$numberofvariables
    } else {
      this$boxnbpointseff <- this$boxnbpoints
    }
    hasbounds <- optimbase.hasbounds(this=this$optbase)
    if (!hasbounds)
      stop('neldermead_updatesimp: Randomized bounds initial simplex is not available without bounds.',
           call.=FALSE)
    simplex0 <- optimsimplex.destroy(this=simplex0)
    tmp <- optimsimplex(method='randbounds',x0=transpose(xopt),fun=costf.transposex,
                            boundsmin=this$optbase$boundsmin,boundsmax=this$optbase$boundsmax,
                            nbve=this$boxnbpointseff,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }

  #
  # Scale the simplex into the bounds and the nonlinear inequality constraints, if any.
  # Caution !
  # The initial simplex may be computed with an 'axis-by-axis' simplex,
  # so that it does not satisfies bounds constraints.
  # The scaling should therefore take into accounts for bounds.
  # TODO : project vertices into bounds
  #
  nbve <- optimsimplex.getnbve(this=simplex0)
  this <- neldermead.log(this=this,msg='Before scaling:')
  str <- optimsimplex.tostring(x=simplex0)
  for (i in 1:nbve){
    this <- neldermead.log(this=this,msg=str[i])
  }
  hasbounds <- optimbase.hasbounds(this=this$optbase)
  hasnlcons <- optimbase.hasnlcons(this=this$optbase)
  if (hasbounds | hasnlcons) {
    this <- neldermead.log(this=this,
                           msg='Scaling initial simplex into nonlinear inequality constraints...')
    nbve <- optimsimplex.getnbve(this=simplex0)
    for (ive in 1:nbve){
      # x is a row vector
      x <- optimsimplex.getx(this=simplex0,ive=ive)
      this <- neldermead.log(this=this,msg=sprintf('Scaling vertex #%d/%d at [%s]... ',ive,nbve,strvec(x)))
      # Transpose x because xopt is a column vector : xp is now a column vector
      tmp <- scaleinconstraints(this=this,x=transpose(x),xref=xopt)
        this <- tmp$this
        status <- tmp$isscaled
        xp <- tmp$p
      rm(tmp) 
      if (!status)
        stop(sprintf('neldermead_updatesimp: Impossible to scale the vertex #%d/%d at [%s] into inequality constraints',
                     ive,nbve,strvec(x)),
             call.=FALSE) 
      if (any(x!=transpose(xp))){
        index <- 2
        tmp <- optimbase.function(this=this$optbase,x=xp,index=index)
          this$optbase <- tmp$this
          fv <- tmp$f
          index <- tmp$index
          if (hasnlcons) c <- tmp$c
        rm(tmp)
        # Transpose xp because optimsimplex takes row coordinate vectors.
        simplex0 <- optimsimplex.setve(this=simplex0,ive=ive,fv=fv,x=transpose(xp))
      }
    }
  }
  this <- neldermead.log(this=this,msg='After scaling:')
  str <- optimsimplex.tostring(x=simplex0)
  for (i in 1:nbve){
    this <- neldermead.log(this=this,msg=str[i])
  }
  this$simplex0 <- optimsimplex.destroy(this=this$simplex0)
  this$simplex0 <- simplex0
  this$simplexsize0 <- optimsimplex.size(this=simplex0)
  this$simplex0 <- optimsimplex.sort(this=this$simplex0)

  return(this)

}

