# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the optimbase component
# originally written by Michael Baudin for Scilab.

optimbase.checkcostfun <- function(this=NULL){

  if (length(this$x0)==0)
    stop('optimbase.checkcostfun: Cannot check cost function when x0 is empty',
         call.=FALSE)
  
  #
  # If there are nonlinear constraints and no derivatives, check that the index is correctly managed.
  #
  if ((this$nbineqconst>0) & (!this$withderivatives)){
    # index: 1
    tmp <- try(optimbase.function(this=this,x=this$x0,index=1))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,1)',
        call.=FALSE)
    } else {
      this <- tmp$this
    }
    # index: 2
    tmp <- try(optimbase.function(this=this,x=this$x0,index=2))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,2)',
        call.=FALSE)
    } else {
      this <- tmp$this
      this <- optimbase.checkshape(this=tmp$this,varname='f',data=tmp$f,
        index=2,expectednrows=1,expectedncols=1)
    }
    # index: 5
    tmp <- try(optimbase.function(this=this,x=this$x0,index=5))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,5)',
        call.=FALSE)
    } else {
      this <- optimbase.checkshape(this=tmp$this,varname='c',data=tmp$c,
        index=5,expectednrows=1,expectedncols=tmp$this$nbineqconst)
    }
    # index: 6
    tmp <- try(optimbase.function(this=this,x=this$x0,index=6))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,6)',
        call.=FALSE)
     } else {
       this <- optimbase.checkshape(this=tmp$this,varname='f',data=tmp$f,
         index=6,expectednrows=1,expectedncols=1)
       this <- optimbase.checkshape(this=this,varname='c',data=tmp$c,
         index=6,expectednrows=1,expectedncols=this$nbineqconst)
     }
  }
  #
  # If there are no nonlinear constraints and no derivatives, check that the index is correctly managed.
  #
  if ((this$nbineqconst==0) & (!this$withderivatives)){
    # index: 1
    tmp <- try(optimbase.function(this=this,x=this$x0,index=1))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,1)',
        call.=FALSE)
    } else {
      this <- tmp$this
    }
    # index: 2
    tmp <- try(optimbase.function(this=this,x=this$x0,index=2))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,2)',
        call.=FALSE)
    } else {
      this <- optimbase.checkshape(this=tmp$this,varname='f',data=tmp$f,
        index=2,expectednrows=1,expectedncols=1)
    }
  }
  #
  # If there are no nonlinear constraints and derivatives, check that the index is correctly managed.
  #
  if ((this$nbineqconst==0) & (this$withderivatives)){
    cmd <- paste('tmp <- optimbase.function(this=this,x=this$x0,index=index)',
                 '  this <- tmp$this',
                 '  f <- tmp$f',
                 '  g <- tmp$g',
                 '  index <- tmp$index',
                 'rm(tmp)',
                 sep='\n')
    # index: 1
    tmp <- try(optimbase.function(this=this,x=this$x0,index=1))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,1)',
        call.=FALSE)
    } else {
      this <- tmp$this
    }
    # index: 2
    tmp <- try(optimbase.function(this=this,x=this$x0,index=2))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,2)',
        call.=FALSE)
    } else {
      this <- optimbase.checkshape(this=tmp$this,varname='f',data=tmp$f,
        index=2,expectednrows=1,expectedncols=1)
    }
    # index: 3
    tmp <- try(optimbase.function(this=this,x=this$x0,index=3))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,3)',
        call.=FALSE)
    } else {
      this <- optimbase.checkshape(this=tmp$this,varname='g',data=tmp$g,
        index=3,expectednrows=1,expectedncols=tmp$this$numberofvariables)
    }
    # index: 4
    tmp <- try(optimbase.function(this=this,x=this$x0,index=4))
    if (class(tmp)=='try-error'){
      stop('optimbase.checkcostfun: Cannot evaluate cost function from costf(x0,4)',
        call.=FALSE)
    } else {
      this <- optimbase.checkshape(this=tmp$this,varname='f',data=tmp$f,
        index=4,expectednrows=1,expectedncols=1)
      this <- optimbase.checkshape(this=this,varname='g',data=tmp$g,
        index=4,expectednrows=1,expectedncols=this$numberofvariables)
    }
  }

  return(this)
}

