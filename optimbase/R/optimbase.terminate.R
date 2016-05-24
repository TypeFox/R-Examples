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

optimbase.terminate <- function(this=NULL,previousfopt=NULL,currentfopt=NULL,
                                previousxopt=NULL,currentxopt=NULL){

  terminate <- FALSE
  status <- 'continue'
  if (this$verbose==TRUE) this <- optimbase.stoplog(this=this,msg=sprintf('  > Termination ?'))

  #
  # Criteria #1 : maximum number of iterations
  #
  if (!terminate){
    if (this$verbose==TRUE)
      this <- optimbase.stoplog(this=this,
                                msg=sprintf('  > iterations=%d >= maxiter=%d',
                                            this$iterations,this$maxiter))
    if (this$iterations>=this$maxiter){
      terminate <- TRUE
      status <- 'maxiter'
    }
  }
  
  #
  # Criteria #2 : maximum number of call to function
  #
  if (!terminate){
    if (this$verbose==TRUE)
      this <- optimbase.stoplog(this=this,
                                msg=sprintf('  > funevals=%d >= maxfunevals=%d',
                                            this$funevals,this$maxfunevals))
    if (this$funevals>=this$maxfunevals){
      terminate <- TRUE
      status <- 'maxfuneval'
    }
  }
  
  #
  # Criteria #3 : tolerance on function
  # Note :
  #   This termination criteria works well in the special case where the function
  #   value at optimum is several order of magnitude smaller
  #   than the initial function value (ie f(x0)).
  #   This is the case when the function value at optimum is zero.
  #   When the function value at optimum is non-zero, or if the
  #   initial function value is strictly positive (e.g. f(x0)=10)
  #   and the optimum function value is strictly negative (e.g. f(x*)=-10),
  #   that criteria fails miserably.
  #
  if (!terminate){
    if (this$tolfunmethod){
      tolfr <- this$tolfunrelative
      tolfa <- this$tolfunabsolute
      acfopt <- abs(currentfopt)
      apfopt <- abs(previousfopt)
      if (this$verbose==TRUE)
        this <- optimbase.stoplog(this=this,
                                  msg=sprintf('  > abs(currentfopt)=%e < tolfunrelative * abs(previousfopt) + tolfunabsolute=%e',
                                              acfopt,tolfr*apfopt+tolfa))
      if (acfopt<tolfr*apfopt+tolfa){
        terminate <- TRUE
        status <- 'tolf'
      }
    }
  }
  
  #
  # Criteria #4 : tolerance on x
  # Note
  # What means a relative error on x ?
  # Notes: if xn and xn+1 are very close to xopt and xopt different from 0,
  # the relative error between xn and xn+1 is small.
  # But if xopt, xn and xn+1 are close to 0, the relative error may be a
  # completely wrong criteria. The absolute tolerance should be used in this case.
  #
  if (!terminate){
    if (this$tolxmethod){
      normdelta <- norm(currentxopt - previousxopt)
      normold <- norm(currentxopt)
      tolxr <- this$tolxrelative
      tolxa <- this$tolxabsolute
      if (this$verbose==TRUE)
        this <- optimbase.stoplog(this=this,
                                  msg=sprintf('  > e(x)=%e < %e * %e + %e',
                                              normdelta,tolxr,normold,tolxa))

      if (normdelta<tolxr*normold+tolxa){
        terminate <- TRUE
        status <- 'tolx'
      }
    }
  }
  if (this$verbose==TRUE)
    this <- optimbase.stoplog(this=this,
                              msg=sprintf('  > Terminate = %s, status = %s',
                                          terminate,status))

  varargout <- list(this=this,terminate=terminate,status=status)

  return(varargout)

}

