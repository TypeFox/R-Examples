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
# This source code is a R port of the optimsimplex component
# originally written by Michael Baudin for Scilab.

optimsimplex.deltafv <- function(this=NULL){

  if (size(this$fv,2)!=1)
    stop(sprintf('optimsimplex.deltafv: this$fv is expected to be a column matrix, but current shape is %d x %d.',
                 size(this$fv,1),size(this$fv,2)),
         call.=FALSE)
  df <- this$fv[2:this$nbve,1] - this$fv[1,1]*ones(this$nbve-1,1)
  return(df)
  
}

