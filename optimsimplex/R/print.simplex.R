# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

print.simplex <- function(x=NULL,...){
  
  cat(sprintf('Dimension: n=%d\n',x$n))
  cat(sprintf('Number of vertices: nbve=%d\n',x$nbve))
  
  str <- optimsimplex.tostring(x)
  
  for (k in 1:x$nbve){
    cat(sprintf('  %s',str[k]))
  }
  
}

