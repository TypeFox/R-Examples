#
# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

fminbnd.function <- function(x=NULL,index=NULL,fmsfundata=NULL){

  fminbnd <- list(f=fmsfundata$Fun(x),
                  index=index,
                  this=list(costfargument=fmsfundata))
  return(fminbnd)
  
}

