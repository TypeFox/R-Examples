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

neldermead.termination <- function(this=NULL,fvinitial=NULL,oldfvmean=NULL,
                                   newfvmean=NULL,previousxopt=NULL,
                                   currentxopt=NULL,simplex=NULL){

  terminate <- FALSE
  status <- 'continue'
  verbose <- optimbase.get(this=this$optbase,key='verbose')
  #
  # Termination Criteria from parent optimization class
  #
  tmp <- optimbase.terminate(this=this$optbase,previousfopt=fvinitial,
                             currentfopt=newfvmean,previousxopt=previousxopt,
                             currentxopt=currentxopt)
    this$optbase <- tmp$this
    terminate <- tmp$terminate
    status <- tmp$status
  rm(tmp)

  #
  # Criteria #6 : simplex absolute + relative size
  #
  if (!terminate){
    if (this$tolsimplexizemethod){
      ssize <- optimsimplex.size(this=simplex,method='sigmaplus')
      tolsa <- this$tolsimplexizeabsolute
      tolsr <- this$tolsimplexizerelative
      ssize0 <- this$simplexsize0
      if (verbose==TRUE)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('  > simplex size=%e < %e + %e * %e',
                                                      ssize,tolsa,tolsr,ssize0))
      if (ssize<(tolsa+tolsr*ssize0)){
        terminate <- TRUE
        status <- 'tolsize'
      }
    }
  }

  #
  # Criteria #7 : simplex absolute size + difference in function values (Matlab-like)
  #
  if (!terminate){
    if (this$tolssizedeltafvmethod){
      ssize <- optimsimplex.size(this=simplex,method='sigmaplus')
      if (verbose==TRUE)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('  > simplex size=%e < %e',
                                                      ssize,this$tolsimplexizeabsolute))
      shiftfv <- abs(optimsimplex.deltafvmax(this=simplex))
      if (verbose==TRUE)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('  > abs(fv(n+1) - fv(1))=%e < toldeltafv=%e',
                                                      shiftfv,this$toldeltafv))
      if ((ssize<this$tolsimplexizeabsolute) & (shiftfv<this$toldeltafv)){
        terminate <- TRUE
        status <- 'tolsizedeltafv'
      }
    }
  }
  
  #
  # Criteria #8 : Kelley stagnation, based on
  # a sufficient decrease condition
  #
  if (!terminate){
    if (this$kelleystagnationflag){
      tmp <- optimsimplex.gradientfv(this=simplex,fun=neldermead.costf,method='forward',data=this)
        sg <- tmp$g
        this <- tmp$data
      rm(tmp)
      nsg <- transpose(sg)%*%sg
      if (verbose==TRUE){
        sgstr <- strvec(sg)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('Test Stagnation: nsg = %e, sg = [%s]',
                                                      nsg,sgstr))
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('Test Stagnation: newfvmean=%e >= oldfvmean=%e - %e * %e',
                                                      newfvmean,oldfvmean,this$kelleyalpha,nsg))
      }
      if (newfvmean>=(oldfvmean-this$kelleyalpha*nsg)){
        terminate <- TRUE
        status <- 'kelleystagnation'
      }
    }
  }

  #
  # Criteria #9 : Box termination criteria
  # The number of consecutive time that an absolute tolerance on
  # function value is met.
  # From Algorithm 454, the tolerance is the difference between the
  # max and the min function values in the simplex.
  #
  if (!terminate){
    if (this$boxtermination){
      shiftfv <- abs(optimsimplex.deltafvmax(this=simplex))
      if (verbose==TRUE)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('Test Box : shiftfv=%e < boxtolf=%e',
                                                      shiftfv,this$boxtolf))
      if (shiftfv<this$boxtolf){
        this$boxkount <- this$boxkount + 1
        if (verbose==TRUE)
          this$optbase <- optimbase.stoplog(this=this$optbase,
                                            msg=sprintf('Test Box : boxkount=%d == boxnbmatch=%d',
                                                        this$boxkount,this$boxnbmatch))
        if (this$boxkount==this$boxnbmatch){
          terminate <- TRUE
          status <- 'tolboxf'
        }
      } else {
        this$boxkount <- 0
      }
    }
  }

  #
  # Criteria #10 : variance of function values
  #
  if (!terminate){
    if (this$tolvarianceflag){
      var <- optimsimplex.fvvariance(this=simplex)
      if (verbose==TRUE)
        this$optbase <- optimbase.stoplog(this=this$optbase,
                                          msg=sprintf('Test tolvariance: %e < %e',
                                                      var,this$tolabsolutevariance))
      if (var<(this$tolrelativevariance*this$variancesimplex0+this$tolabsolutevariance)){
        terminate <- TRUE
        status <- 'tolvariance'
      }
    }
  }

  #
  # Criteria #11 : user-defined criteria
  #
  if (!terminate){
    if (this$myterminateflag){
      tmp <- this$myterminate(this,simplex)   # this$myterminate is supposedly a function, but this functionnality does not seems to be used in the current neldermead algorithm
        this <- tmp$this
        term <- tmp$term
        stat <- tmp$stat
      rm(tmp)
      if (term){
        terminate <- term
        status <- stat
      }
    }
  }
  if (verbose==TRUE)
    this$optbase <- optimbase.stoplog(this=this$optbase,
                                      msg=sprintf('  > Terminate = %s, status = %s',
                                                  terminate,status))

  varargout <- list(this=this,terminate=terminate,status=status)

  return(varargout)

}

