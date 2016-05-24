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

optimbase.get <- function(this=NULL,key=NULL){

  if (!any(key==c('verbose','x0','fx0','xopt','fopt','tolfunabsolute',
      'tolfunrelative','tolfunmethod','tolxabsolute','tolxrelative',
      'tolxmethod','maxfunevals','funevals','maxiter','iterations','function',
      'status','historyfopt','historyxopt','verbosetermination','outputcommand',
      'outputcommandarg','numberofvariables','storehistory','costfargument',
      'boundsmin','boundsmax','nbineqconst','logfile','logfilehandle',
      'logstartup','withderivatives'))){
    stop(sprintf('optimbase.get: Unknown key %s',key),
      call.=FALSE)
  }
  
  if ((key=='historyxopt' | key=='historyfopt') & (!this$storehistory)){
    stop('optimbase.get: History disabled ; enable storehistory option.',
      call.=FALSE)
  }
  
  return(this[[key]])
  
}

