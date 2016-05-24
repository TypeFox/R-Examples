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

optimbase.checkshape <- function(this=NULL,varname=NULL,data=NULL,
                                 index=NULL,expectednrows=NULL,expectedncols=NULL){

  if (size(data,1)!=expectednrows)
    stop(sprintf('optimbase.checkshape: The matrix %s from costf(x0,%d) has %d rows, instead of %d.',
                  varname,index,size(data,1),expectednrows,),
          call.=FALSE)

  if (size(data,2)!=expectedncols)
    stop(sprintf('optimbase.checkshape: The matrix %s from costf(x0,%d) has %d columns, instead of %d.',
                  varname,index,size(data,2),expectedncols),
          call.=FALSE)
  return(this)
}

