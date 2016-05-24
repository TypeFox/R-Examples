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

optimbase.outstruct <- function(this=NULL){

  if (is.null(this$outputcommand)){
    stop('optimbase.outstruct: No output command is defined.',call.=FALSE)
  } else {
    data <- list(x=this$xopt,
                 fval=this$fopt,
                 iteration=this$iterations,
                 funccount=this$funevals)
    class(data) <- 'optimbase.outputdata'
  }
  
  return(data)
  
}

