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

neldermead.startup <- function(this=NULL){

  # 0. Check that the cost function is correctly connected
  # Note: this call to the cost function is not used, but helps the
  # user while he is tuning his object.
  if (this$checkcostfunction)
    this$optbase <- optimbase.checkcostfun(this=this$optbase)

  # 1. If the problem has bounds, check that they are consistent
  hasbounds <- optimbase.hasbounds(this=this$optbase)
  if (hasbounds){
    tmp <- optimbase.checkbounds(this=this$optbase)
    if (!tmp$isok)
      stop(sprintf('neldermead.startup: %s',tmp$errmsg),
           call.=FALSE)
    rm(tmp)
  }

  # 2. Get the initial guess and compute the initial simplex
  x0 <- optimbase.get(this=this$optbase,key='x0')
  if (!any(this$simplex0method==c('given','axes','spendley','pfeffer','randbounds')))
    stop(sprintf('neldermead.startup: Unknown value %s for -simplex0method option.',
                 this$simplex0method),
         call.=FALSE)

  if (this$simplex0method=='given'){
    tmp <- optimsimplex(coords=this$coords0,fun=costf.transposex,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$simplex0method=='axes'){
    tmp <- optimsimplex(method='axes',x0=transpose(x0),fun=costf.transposex,
                            len=this$simplex0length,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$simplex0method=='spendley'){
    tmp <- optimsimplex(method='spendley',x0=transpose(x0),fun=costf.transposex,
                            len=this$simplex0length,data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$simplex0method=='pfeffer'){
    tmp <- optimsimplex(method='pfeffer',x0=transpose(x0),fun=costf.transposex,
                            deltausual=this$simplex0deltausual,deltazero=this$simplex0deltazero,
                            data=this)
      simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }
  if (this$simplex0method=='randbounds'){
    if (this$boxnbpoints=='2n'){
      this$boxnbpointseff <- 2 * this$optbase$numberofvariables
    } else {
      this$boxnbpointseff <- this$boxnbpoints
    }
    if (!hasbounds)
      stop('neldermead.startup: Randomized bounds initial simplex is not available without bounds.',
           call.=FALSE)
    tmp <- optimsimplex(method='randbounds',x0=transpose(x0),fun=costf.transposex,
                            boundsmin=this$optbase$boundsmin,boundsmax=this$optbase$boundsmax,
                            nbve=this$boxnbpointseff,data=this)
       simplex0 <- tmp$newobj
      this <- tmp$data
    rm(tmp)
  }

  #
  # 3. Scale the initial simplex into the bounds and the nonlinear inequality constraints, if any
  #
  hasnlcons <- optimbase.hasnlcons(this=this$optbase)
  if (hasbounds | hasnlcons){
    # Check that initial guess is feasible
    tmp <- optimbase.isfeasible(this=this$optbase,x=x0)
      this$optbase <- tmp$this
      isfeasible <- tmp$isfeasible
    rm(tmp)
    if (isfeasible!=1)
      stop(sprintf('neldermead.startup: Initial guess (%s) is not feasible.',strvec(x0)),
           call.=FALSE)
    this <- neldermead.log(this=this,
                           msg='Scaling initial simplex into nonlinear inequality constraints...')

    if (!any(this$scalingsimplex0==c('tox0','tocenter')))
      stop(sprintf('neldermead.startup: Unknown value %s for -scalingsimplex0 option',
                   this$scalingsimplex0),
           call.=FALSE)
    if (this$scalingsimplex0=='tox0')
      simplex0 <- neldermead.scaletox0(this=this,simplex0=simplex0)
    if (this$scalingsimplex0=='tocenter')
      simplex0 <- neldermead.scaletocenter(this=this,simplex0=simplex0)
  }

  #
  # 4. Store the simplex
  #
  this$simplex0 <- optimsimplex.destroy(this=this$simplex0)
  this$simplex0 <- simplex0
  this$simplexsize0 <- optimsimplex.size(this=simplex0)

  # 5. Store initial data into the base optimization component
  fx0 <- optimsimplex.getfv (this=this$simplex0,ive=1)
  this$optbase <- optimbase.set(this=this$optbase,key='fx0',value=fx0)
  this$optbase <- optimbase.set(this=this$optbase,key='xopt',value=x0)
  this$optbase <- optimbase.set(this=this$optbase,key='fopt',value=fx0)
  this$optbase <- optimbase.set(this=this$optbase,key='iterations',value=0)

  # 6. Initialize the termination criteria
  this <- neldermead.termstartup(this=this)

  return(this)

}

