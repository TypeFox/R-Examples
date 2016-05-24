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

optimbase.histset <- function(this=NULL,iter=NULL,key=NULL,value=NULL){

  if (!this$storehistory)
    stop('optimbase.histset: History disabled ; turn on storehistory option.',
         call.=FALSE)

  if (iter<1)
    stop('optimbase.histset: Negative iteration index are not allowed.',
         call.=FALSE)

  if (!any(key==c('historyxopt','historyfopt')))
    stop(sprintf('optimbase.histset: Unknown key %s',key),
         call.=FALSE)

  if (key=='historyxopt') this$historyxopt[[iter]] <- value
  if (key=='historyfopt') this$historyfopt[iter] <- value

  return(this)
  
}

