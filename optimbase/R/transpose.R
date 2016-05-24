# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

transpose <- function(object=NULL){
  if (is.null(object))
    stop('transpose: input is not a matrix or a vector.',
         call.=FALSE)
  if (!is.matrix(object)){
    if (is.atomic(object)){
      mat <- matrix(object,ncol=length(object),byrow=TRUE)
    } else {
      stop('transpose: input is not a matrix or a vector.',
           call.=FALSE)
    }
  } else {
    mat <- object
  }
  return(t(mat))
}

