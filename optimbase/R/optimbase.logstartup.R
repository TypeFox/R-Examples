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

optimbase.logstartup <- function(this=NULL){

  if (this$logstartup){
    stop('optimbase.logstartup: Logging already started.',
         call.=FALSE)
  } else {
    this$logstartup <- TRUE
    if (this$logfile != ''){
      if (this$logfilehandle != 0){
        stop('optimbase.logstartup: Log file handle non zero while starting up the logging.',
             call.=FALSE)
      }
      this$logfilehandle <- 1
      write(paste('Optimbase ready for logging at',Sys.time()),
            file=this$logfile,
            append=TRUE)
    }
  }
  
  return(this)
  
}

