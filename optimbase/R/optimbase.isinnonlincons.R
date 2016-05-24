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
# This source code is a R port of the optimbase component
# originally written by Michael Baudin for Scilab.

optimbase.isinnonlincons <- function(this=NULL,x=NULL){

  isfeasible <- TRUE
  if (this$nbineqconst > 0){
    if (this$withderivatives){
      tmp <- optimbase.function(this=this,x=x,index=2)
        this <- tmp$this
        f <- tmp$f
        g <- tmp$g
        c <- tmp$c
        gc <- tmp$gc
        index <- tmp$index
      rm(tmp)
    } else {
      tmp <- optimbase.function(this=this,x=x,index=5)
        this <- tmp$this
        f <- tmp$f
        c <- tmp$c
        index <- tmp$index
      rm(tmp)
    }
    for (ic in 1:this$nbineqconst){
      if (c[ic] < 0.0){
        this <- optimbase.log(this=this,
                              msg=sprintf('Inequality constraint #%d/%d is not satisfied for x',
                                          ic,this$nbineqconst))
        isfeasible <- FALSE
        break
      }
    }
  }

  varargout <- list(this=this,isfeasible=isfeasible)

  return(varargout)

}

