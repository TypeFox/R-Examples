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

optimbase.function <- function(this=NULL,x=NULL,index=NULL){

  if (is.null(this) | is.null(x) | is.null(index))
    stop('optimbase.function: one or more input argument is missing.',
         call.=FALSE)

  if (is.null(this$fun))
    stop('optimbase.function: Empty function (use -function option).',
         call.=FALSE)

  #Start evaluating the cost function
  this$funevals <- this$funevals + 1;

  # Setup output
  varargout <- list(this=this,f=NULL,g=NULL,c=NULL,gc=NULL,index=index)

  if (this$verbose==TRUE)
    this <- optimbase.log(this=this,
                          msg=sprintf('Function Evaluation #%d at [%s]',
                                      this$funevals,
                                      paste(x,collpase=' ')))

  if (this$withderivatives){
    if (this$nbineqconst==0){
      if (is.null(this$costfargument) ||
        class(this$costfargument)!='optimbase.functionargs'){
        tmp <- this$fun(x=x,index=index)
          varargout$f <- tmp$f
          varargout$g <- tmp$g
          varargout$index <- tmp$index
        rm(tmp)
      } else {
        tmp <- this$fun(x=x,index=index,fmsfundata=this$costfargument)
          varargout$f <- tmp$f
          varargout$g <- tmp$g
          varargout$index <- tmp$index
          varargout$this$costfargument <- tmp$this$costfargument
        rm(tmp)
      }
    } else {
      if (is.null(this$costfargument) ||
        class(this$costfargument)!='optimbase.functionargs'){
        tmp <- this$fun(x=x,index=index)
          varargout$f <- tmp$f
          varargout$g <- tmp$g
          varargout$c <- tmp$c
          varargout$gc <- tmp$gc
          varargout$index <- tmp$index
        rm(tmp)
      } else {
        tmp <- this$fun(x=x,index=index,fmsfundata=this$costfargument)
          varargout$f <- tmp$f
          varargout$g <- tmp$g
          varargout$c <- tmp$c
          varargout$gc <- tmp$gc
          varargout$index <- tmp$index
          varargout$this$costfargument <- tmp$this$costfargument
        rm(tmp)
      }
    }
  } else {
    if (this$nbineqconst==0){
      if (is.null(this$costfargument) ||
        class(this$costfargument)!='optimbase.functionargs'){
        tmp <- this$fun(x=x,index=index)
          varargout$f <- tmp$f
          varargout$index <- tmp$index
        rm(tmp)
      } else {
        tmp <- this$fun(x=x,index=index,fmsfundata=this$costfargument)
          varargout$f <- tmp$f
          varargout$index <- tmp$index
          varargout$this$costfargument <- tmp$this$costfargument
        rm(tmp)
      }
    } else {
      if (is.null(this$costfargument) ||
        class(this$costfargument)!='optimbase.functionargs'){
        tmp <- this$fun(x,index)
          varargout$f <- tmp$f
          varargout$c <- tmp$c
          varargout$index <- tmp$index
        rm(tmp)
      } else {
        tmp <- this$fun(x=x,index=index,fmsfundata=this$costfargument)
          varargout$f <- tmp$f
          varargout$c <- tmp$c
          varargout$index <- tmp$index
          varargout$this$costfargument <- tmp$this$costfargument
        rm(tmp)
      }
    }
  }

  return(varargout)

}

