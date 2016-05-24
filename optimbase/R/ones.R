# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

ones <- function(nx = 1, ny = nx){

  if (is.null(nx) | is.null(ny))
    stop('ones: nx or ny is NULL',
         call.=FALSE)

  if (!is.numeric(nx) | !is.numeric(ny))
    stop('ones: nx or ny is not numeric',
         call.=FALSE)

  if (length(nx)!=1 | length(ny)!=1)
    stop('ones: nx or ny is not of length 1',
         call.=FALSE)

  if (nx<=0 | ny<=0)
    stop('ones: nx or ny is less or equal to 0',
         call.=FALSE)

  nx <- as.integer(nx)
  ny <- as.integer(ny)

  return(matrix(1, nrow=nx, ncol=ny))
  
}


zeros <- function(nx = 1, ny = nx){

  if (is.null(nx) | is.null(ny))
    stop('zeros: nx or ny is NULL',
         call.=FALSE)

  if (!is.numeric(nx) | !is.numeric(ny))
    stop('zeros: nx or ny is not numeric',
         call.=FALSE)

  if (length(nx)!=1 | length(ny)!=1)
    stop('zeros: nx or ny is not of length 1',
         call.=FALSE)

  if (nx<=0 | ny<=0)
    stop('zeros: nx or ny is less or equal to 0',
         call.=FALSE)

  nx <- as.integer(nx)
  ny <- as.integer(ny)

  return(matrix(0, nrow=nx, ncol=ny))

}