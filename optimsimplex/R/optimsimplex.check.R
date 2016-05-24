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

optimsimplex.check <- function(this=NULL){

  nx1 <- size(this$x,1)
  nx2 <- size(this$x,2)
  if (this$nbve!=0 & nx1!=this$nbve)
    stop(sprintf('optimsimplex.check: Number of rows of x is %d, which is different from number of vertices = %d.',
                 nx1,this$nbve),
         call.=FALSE)
  if (this$n!=0 & nx2!=this$n)
    stop(sprintf('optimsimplex.check: Number of columns of x is %d, which is different from dimension = %d.',
                 nx2,this$n),
         call.=FALSE)

  nf1 <- size(this$fv,1)
  nf2 <- size(this$fv,2)
  if (this$n!=0 & nf1!=this$nbve)
    stop(sprintf('optimsimplex.check: Number of rows of fv is %d, which is different from number of vertices = %d.',
                 nf1,this$nbve),
         call.=FALSE)
  if (this$nbve!=0 & nf2!=1)
    stop(sprintf('optimsimplex.check: Number of columns of fv is %d, which is different from 1.',
                 nf2),
         call.=FALSE)

}

