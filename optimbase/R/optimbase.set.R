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

optimbase.set <- function(this=NULL,key=NULL,value=NULL){
  
  if (!any(key==c('verbose','x0','fx0','xopt','fopt','tolfunabsolute',
      'tolfunrelative','tolfunmethod','tolxabsolute','tolxrelative',
      'tolxmethod','maxfunevals','funevals','maxiter','iterations','function',
      'status','historyfopt','historyxopt','verbosetermination','outputcommand',
      'outputcommandarg','numberofvariables','storehistory','costfargument',
      'boundsmin','boundsmax','nbineqconst','logfile','logfilehandle',
      'logstartup','withderivatives'))){
    stop(sprintf('optimbase.set: Unknown key %s',key),
      call.=FALSE)
  }
  
  if (key=='iterations'){
    this$iterations <- value
  }
  else if (key=='xopt'){
    this$xopt <- value
  }
  else if (key=='fopt'){
    this$fopt <- value
  }
  else if (key=='historyxopt'){
    if (!this$storehistory){
      stop('optimbase.set: History disabled ; turn on -storehistory option.',
        call.=FALSE)
    } else {
      this$historyxopt <- value
    }
  }
  else if (key=='historyfopt'){
    if (!this$storehistory){
      stop('optimbase.set: History disabled ; turn on -storehistory option.',
        call.=FALSE)
    } else {
      this$historyfopt <- value
    }
  }
  else if (key=='fx0'){
    this$fx0 <- value
  }
  else if (key=='status'){
    this$status <- value
  }
  else if (key=='verbose') {
    assert.classboolean(var=value,varname='value',ivar=3)
    this$verbose <- value
  }
  else if (key=='verbosetermination') {
    assert.classboolean(var=value,varname='value',ivar=3)
    this$verbosetermination <- value
  }
  else if (key=='logfile') {
    if (this$logstartup) this <- optimbase.logshutdown(this=this)
    this$logfile <- value
    this <- optimbase.logstartup(this=this)
  }
  else if (key=='x0') {
    assert.classreal(var=value,varname='value',ivar=3)
    if (size(value,2)!=1)
      stop(sprintf('optimbase.set: The x0 vector is expected to be a column matrix, but current shape is %d x %d.',
          size(value,1),size(value,2)),
        call.=FALSE)
    this$x0 <- value
  }
  else if (key=='maxfunevals') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$maxfunevals <- value
  }
  else if (key=='funevals'){
    assert.classreal(var=value,varname='value',ivar=3)
    this$funevals <- value
  }
  else if (key=='maxiter') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$maxiter <- value
  }
  else if (key=='tolfunabsolute') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$tolfunabsolute <- value
  }
  else if (key=='tolfunrelative') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$tolfunrelative <- value
  }
  else if (key=='tolxabsolute') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$tolxabsolute <- value
  }
  else if (key=='tolxrelative') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$tolxrelative <- value
  }
  else if (key=='tolxmethod') {
    assert.classboolean(var=value,varname='value',ivar=3)
    if(!is.logical(value)){
      unknownValueForOption(value=value,optionname='tolxmethod')
    } else {
      this$tolxmethod <- value
    }
  }
  else if (key=='tolfunmethod') {
    assert.classboolean(var=value,varname='value',ivar=3)
    if(!is.logical(value)){
      unknownValueForOption(value=value,optionname='tolfunmethod')
    } else {
      this$tolfunmethod <- value
    }
  }
  else if (key=='function') {
    assert.classfunction(var=value,varname='value',ivar=3)
    this$fun <- value
  }
  else if (key=='outputcommand') {
    assert.classfunction(var=value,varname='value',ivar=3)
    this$outputcommand <- value
  }
  else if (key=='outputcommandarg') {
    if (!is(value)=='optimbase.outputargs'){
      stop(paste('optimbase.set: the outputcommandarg argument must be an',
          'optimbase.outputargs object.'),call.=FALSE)
    } else {
      this$outputcommandarg <- value
    }
  }
  else if (key=='numberofvariables') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$numberofvariables <- value
  }
  else if (key=='storehistory') {
    assert.classboolean(var=value,varname='value',ivar=3)
    this$storehistory <- value
  }
  else if (key=='costfargument') {
    if (!is(value)=='optimbase.functionargs'){
      stop(paste('optimbase.set: the costargument argument must be an',
          'optimbase.functionargs object.'),call.=FALSE)
    } else {
      this$costfargument <- value
    }
  }
  else if (key=='boundsmin') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$boundsmin <- value
  }
  else if (key=='boundsmax') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$boundsmax <- value
  }
  else if (key=='nbineqconst') {
    assert.classreal(var=value,varname='value',ivar=3)
    this$nbineqconst <- value
  }
  else if (key=='withderivatives') {
    assert.classboolean(var=value,varname='value',ivar=3);
    this$withderivatives <- value
  }
  
  return(this)
  
}

