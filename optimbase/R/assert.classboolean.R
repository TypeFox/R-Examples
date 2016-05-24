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

assert.classboolean <- function(var=NULL,varname=NULL,ivar=NULL){
  if (is.null(var))
    stop('assert.classboolean: no input specified, the variable cannot be checked.',
         call.=FALSE)
  # Generates an error if the given variable is not of type boolean
  if (!is.logical(var)){
    if (is.null(varname) | is.null(ivar)){
      stop(sprintf('assert.classboolean: Expected real variable but got %s instead.',
                   class(var)),
           call.=FALSE)
    } else {
      stop(sprintf('assert.classboolean: Expected logical variable for variable %s at input #%d, but got %s instead.',
                   varname,ivar,class(var)),
           call.=FALSE)
    }
  }
}

