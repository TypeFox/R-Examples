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

optimsimplex.xbar <- function(this=NULL,iexcl=NULL){

  if (is.null(iexcl)) iexcl <- this$nbve

  if (size(iexcl,1)!=1)
    stop(sprintf('optimsimplex.xbar: The exclusion index vector has %d rows instead of 1.',size(iexcl,1)),
         call.=FALSE)

  # Vectorize by making the sum of all, substracting only one vector
  cen <- apply(this$x[1:this$nbve,1:this$n,drop=FALSE],2,sum)
  if (size(this$x[iexcl,1:this$n,drop=FALSE],1)==1){
    tmp <- matrix(this$x[iexcl,1:this$n],ncol=length(1:this$n))
  } else {
    tmp <- this$x[iexcl,1:this$n,drop=FALSE]
  }
  cen <- cen - apply(tmp,2,sum)
  nexcl <- size(iexcl,2)
  cen <- cen/(this$nbve-nexcl)

  return(cen)

}

