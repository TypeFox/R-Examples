# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

size <- function(x=NULL,n=NULL){
  if (is.list(x) & !is.data.frame(x)){
    # x is a list
    return(lapply(x,size))
  } else {
    if (is.null(x)){
      # x is NULL
      if (is.null(n)){
        return(c(0,0))
      } else {
        return(c(0,0)[n])
      }
    } else {
      if (is.null(dim(x))){
        # x is a vector
        if (is.null(n)){
          return(c(1,length(x)))
        } else {
          return(c(1,length(x))[n])
        }
      } else {
        # x is a matrix, an array or a data.frame
        if (is.null(n)){
          return(dim(x))
        } else {
          return(dim(x)[n])
        }
      }
    }
  }
}

