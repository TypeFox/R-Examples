# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of x distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

print.optimbase <- function(x,verbose=FALSE,...){
  
  summary(x,showhistory=verbose)
  
  cat('\nOptimbase Object Definition:\n')
  x[c('numberofvariables','fun','costfargument','x0','fx0','xopt','fopt',
      'boundsmin','boundsmax','funevals','maxfunevals','maxiter','nbineqconst',
      'iterations','status')] <- NULL
  if (verbose) {
    x$historyfopt <- NULL
    x$historyxopt <- NULL
    cat(sprintf('- Optimization Method Based on Derivatives: %s\n',x$withderivatives))
    cat(sprintf('- Verbose logging: %s\n',x$verbose))
    cat(sprintf('- Termination Method on function value: %s\n',x$tolfunmethod))
    cat(sprintf('- Termination Absolute Tolerance on function value: %e\n',x$tolfunabsolute))
    cat(sprintf('- Termination Relative Tolerance on function value: %e\n',x$tolfunrelative))
    cat(sprintf('- Termination Method on x: %s\n',x$tolxmethod))
    cat(sprintf('- Termination Absolute Tolerance on x: %e\n',x$tolxabsolute))
    cat(sprintf('- Termination Relative Tolerance on x: %e\n',x$tolxrelative))
    cat(sprintf('- Verbose Termination: %s\n',x$verbosetermination))
    cat(sprintf('- Verbose Log File: %s\n',x$logfile))
    cat(sprintf('- Verbose Log File Startup Up: %s\n',x$logstartup))
    cat(sprintf('- Store History: %s\n',x$storehistory))
    cat('- Output Command Function:\n')
    print(x$outputcommand)
    cat('- Output Command Function Argument(s):\n')
    print(x$outputcommandarg)
  } else {
    str(x)
  }

}
