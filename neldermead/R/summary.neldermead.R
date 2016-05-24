# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

summary.neldermead <- function(object,showhistory,...){
  
  if(missing('showhistory'))
    showhistory <- FALSE
  
  summary(object$optbase,showhistory=showhistory)
  
  #
  # Simplex info 
  #
  cat('Simplex Information:\n')
  cat('\n- Simplex at Initial Point:\n')
  print(object$simplex0)
  cat('\n- Simplex at Optimal Point:\n')
  print(object$simplexopt)
  
  if (showhistory==TRUE && length(object$historysimplex)>1){
    cat('\n- Simplex History:\n')
    for (iter in 1:length(object$historysimplex)){
      cat(sprintf('* Iteration %d:\n',iter))
      print(object$historysimplex[[iter]])
    }
  }
}