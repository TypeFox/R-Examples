# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

print.vertex <- function(x=NULL,...){
  
  cat(sprintf('Dimension: n=%d\n',x$n))
  ss <- sprintf('%e',x$x[1,1])
  if (x$n>=2){
    for (i in 2:x$n){
      ss <- paste(ss,sprintf('%e',x$x[1,i]))
    }
  }
  cat(sprintf('Vertex: fv=%e, x=%s\n',x$fv,ss))

}

