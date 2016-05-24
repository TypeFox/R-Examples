# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

fminbnd.outputfun <- function(state=NULL,data=NULL,fmsdata=NULL){

  #
  # Compute procedure
  #
  if (!any(data$step==c('init','done','reflection','expansion','insidecontraction',
                        'outsidecontraction','shrink','boxreflection')))
    stop(sprintf('fminbnd.outputfun: Unknown step %s',data$step),
         call.=FALSE)

  if (data$step=='init'){
    if (data$iteration == 0){
      procedure <- ''
    } else {
      procedure <- 'initial simplex'
    }
  }
  if (data$step=='done') procedure <- ''
#  if (data$step=='reflection') procedure <- 'reflect'
#  if (data$step=='expansion') procedure <- 'expand'
#  if (data$step=='insidecontraction') procedure <- 'contract inside'
#  if (data$step=='outsidecontraction') procedure <- 'contract outside'
#  if (data$step=='shrink') procedure <- 'shrink'
  if (data$step=='boxreflection') procedure <- 'reflect (Box)'

  #
  # Display a message
  #
  if (fmsdata$Display=='iter'){
    if (data$fval>=1 & data$fval<1e+8){
      fval <- as.character(signif(data$fval,8))
    }
    if (data$fval>1e-5 & data$fval<1){
      fval <- as.character(round(data$fval,7))
    }
    if (data$fval<=1e-5 | data$fval>=1e+8) {
      fval <- sprintf('%.3e',data$fval)
    }
    if (data$iteration==0){
      data$funccount <- 1
    }
    if (data$step!='done'){
      cat(sprintf('%6s        %5s     %12s         %-20s\n',
                  data$iteration,data$funccount,fval,procedure))
    } else {
      cat('\n') #Why last iteration is not done?
    }
  }
  
  #
  # Process output functions
  #
  optimValues <- list(funccount=data$funccount,
                      fval=data$fval,
                      iteration=data$iteration,
                      procedure=procedure)
  if (length(fmsdata$OutputFcn)!=0){
    if (!any(mode(fmsdata$OutputFcn)==c('function','list')))
      stop('fminsearch: The value of the \'OutputFcn\' option is neither a function nor a list.',
           call.=FALSE)
    if (mode(fmsdata$OutputFcn)=='function'){
      # The output function is a macro
      fmsdata$OutputFcn(data$x,optimValues,state)
    } else {
      # The output function is a list of macros
      for (i in 1:length(fmsdata$OutputFcn)){
        if (!is.list(fmsdata$OutputFcn[i]))
          stop(sprintf(paste('fminsearch: The value of the %dth element of the \'OutputFcn\' option',
                             'is neither a function nor a list.'),i),
               call.=FALSE)
        fmsdata$OutputFcn[i](data$x,optimValues,state)
      }
    }
  }

  # Process plot functions
  if (length(fmsdata$PlotFcns)!=0){
    if (!any(mode(fmsdata$PlotFcns)==c('function','list')))
      stop('fminsearch: The value of the \'PlotFcns\' option is neither a function nor a list.',
           call.=FALSE)
    if (mode(fmsdata$PlotFcns)=='function'){
      # The output function is a macro
      fmsdata$PlotFcns(data$x,optimValues,state)
    } else {
      # The output function is a list of macros
      for (i in 1:length(fmsdata$PlotFcns)){
        if (!is.list(fmsdata$PlotFcns[i]))
          stop(sprintf(paste('fminsearch: The value of the %dth element of the \'PlotFcns\' option',
                             'is neither a function nor a list.'),i),
               call.=FALSE)
        fmsdata$PlotFcns[i](data$x,optimValues,state)
      }
    }
  }
}

