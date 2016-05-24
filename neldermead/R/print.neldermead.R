# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

print.neldermead <- function(x,verbose=FALSE,...){
  
  summary(x,showhistory=verbose)
  
  cat('\nNelder-Mead Object Definition:\n')
  x[c('optbase','simplex0','simplexopt')] <- NULL
  if (verbose) {
    x$historysimplex <- NULL
    cat(sprintf('- Simplex0 Method: %s\n', x$simplex0method)) ; x$simplex0method <-NULL
    cat(sprintf('- Simplex0 Length: %d\n', x$simplex0length)) ; x$simplex0length <-NULL
    cat(sprintf('- Simplex0 Coordinates:\n'))
    print(x$coords0) ; x$coords0 <-NULL
    cat(sprintf('- Simplex0, Pfeffer, Delta Usual: %e\n', x$simplex0deltausual)) ; x$simplex0deltausual <-NULL
    cat(sprintf('- Simplex0, Pfeffer, Delta Zero: %e\n', x$simplex0deltazero)) ; x$simplex0deltazero <-NULL
    cat(sprintf('- Termination Parameters\n'))
    cat(sprintf('  * Termination on simplex size: %s\n', x$tolsimplexizemethod)) ; x$tolsimplexizemethod <-NULL
    cat(sprintf('  * Termination on simplex size, Absolute Tolerance: %e\n', x$tolsimplexizeabsolute)) ; x$tolsimplexizeabsolute <-NULL
    cat(sprintf('  * Termination on simplex size, Relative Tolerance: %e\n', x$tolsimplexizerelative)) ; x$tolsimplexizerelative <-NULL
    cat(sprintf('  * Termination on simplex size, Initial Simplex Size: %d\n', x$simplexsize0)) ; x$simplexsize0 <-NULL
    cat(sprintf('  * Termination on simplex size + Delta of function value: %s\n', x$tolssizedeltafvmethod)) ; x$tolssizedeltafvmethod <-NULL
    cat(sprintf('  * Termination on simplex size + Delta of function value, Absolute Tolerance on Delta F: %e\n', x$toldeltafv)) ; x$toldeltafv <-NULL
    cat(sprintf('  * Termination on Kelley\'s Stagnation: %s\n', x$kelleystagnationflag)) ; x$kelleystagnationflag <-NULL
    cat(sprintf('  * Termination on Kelley\'s Stagnation, Normalization: %s\n', x$kelleynormalizationflag)) ; x$kelleynormalizationflag <-NULL
    cat(sprintf('  * Termination on Kelley\'s Stagnation, Alpha0: %e\n', x$kelleystagnationalpha0)) ; x$kelleystagnationalpha0 <-NULL
    cat(sprintf('  * Termination on Kelley\'s Stagnation, Alpha: %e\n', x$kelleyalpha)) ; x$kelleyalpha <-NULL
    cat(sprintf('  * Termination by Box: %s\n', x$boxtermination)) ; x$boxtermination <-NULL
    cat(sprintf('  * Termination by Box, Absolute Tolerance on Function: %e\n', x$boxtolf)) ; x$boxtolf <-NULL
    cat(sprintf('  * Termination by Box, Maximum Number of Consecutive Match: %d\n', x$boxnbmatch)) ; x$boxnbmatch <-NULL
    cat(sprintf('  * Termination by Box, Current Number of Consecutive Match: %d\n', x$boxkount)) ; x$boxkount <-NULL
    cat(sprintf('  * Termination on Variance: %s\n', x$tolvarianceflag)) ; x$tolvarianceflag <-NULL
    cat(sprintf('  * Termination on Variance, Absolute Tolerance: %e\n', x$tolabsolutevariance)) ; x$tolabsolutevariance <-NULL
    cat(sprintf('  * Termination on Variance, Relative Tolerance: %e\n', x$tolrelativevariance)) ; x$tolrelativevariance <-NULL
    cat(sprintf('  * Termination on Variance, Variance of Initial Simplex: %e\n', x$variancesimplex0)) ; x$variancesimplex0 <-NULL
    cat(sprintf('  * Exit Flag: %s\n', x$exitflag)) ; x$exitflag <-NULL
    # The 2 following elements are currently not used 
    #cat(sprintf('  * Termination, Tolerance on Function Standard Deviation Method: %s\n', x$tolfstdeviationmethod)) ; x$tolfstdeviationmethod <-NULL
    #cat(sprintf('  * Termination, Tolerance on Function Standard Deviation: %s\n', x$tolfstdeviation)) ; x$tolfstdeviation <-NULL
    cat(sprintf('- Algorithms parameters\n'))
    cat(sprintf('  * Method: %s\n', x$method)) ; x$method <-NULL
    cat(sprintf('  * Reflection Factor (rho): %e\n', x$rho)) ; x$rho <-NULL
    cat(sprintf('  * Expansion Factor (chi): %e\n', x$chi)) ; x$chi <-NULL
    cat(sprintf('  * Contraction Factor (gamma): %e\n', x$gamma)) ; x$gamma <-NULL
    cat(sprintf('  * Shrinkage Factor (sigma): %e\n', x$sigma)) ; x$sigma <-NULL
    cat(sprintf('  * Kelley Stagnation: %s\n', x$kelleystagnationflag)) ; x$kelleystagnationflag <-NULL
    cat(sprintf('  * Restart Epsilon: %e\n', x$restarteps)) ; x$restarteps <-NULL
    cat(sprintf('  * Restart Step: %d\n', x$restartstep)) ; x$restartstep <-NULL
    cat(sprintf('  * Restart Maximum: %d\n', x$restartmax)) ; x$restartmax <-NULL
    cat(sprintf('  * Restart Simplex Method: %s\n', x$restartsimplexmethod)) ; x$restartsimplexmethod <-NULL
    cat(sprintf('  * Restart Flag: %s\n', x$restartflag)) ; x$restartflag <-NULL
    cat(sprintf('  * Restart Number: %d\n', x$restartnb)) ; x$restartnb <-NULL
    cat(sprintf('  * Restart Detection Method: %s\n', x$restartdetection)) ; x$restartdetection <-NULL
    cat(sprintf('  * Startup Flag: %s\n', x$startupflag)) ; x$startupflag <-NULL
    cat(sprintf('  * Automatic Checking of Cost Function: %s\n', x$checkcostfunction)) ; x$checkcostfunction <-NULL
    cat(sprintf('  * Box, Number of Points: %s\n', x$boxnbpoints)) ; x$boxnbpoints <-NULL
    cat(sprintf('  * Box, Current Number of Points: %d\n', x$boxnbpointseff)) ; x$boxnbpointseff <-NULL
    cat(sprintf('  * Box, Scaling, Factor: %e\n', x$boxineqscaling)) ; x$boxineqscaling <-NULL
    cat(sprintf('  * Box, Scaling, Method: %s\n', x$scalingsimplex0)) ; x$scalingsimplex0 <-NULL
    cat(sprintf('  * Box, Scaling, Minimum: %e\n', x$guinalphamin)) ; x$guinalphamin <-NULL
    cat(sprintf('  * Box, Bounds Parameter: %e\n', x$boxboundsalpha)) ; x$boxboundsalpha <-NULL
    cat(sprintf('  * Box, Reflection Coefficient: %e\n', x$boxreflect)) ; x$boxreflect <-NULL
    cat(sprintf('  * Greedy: %s\n', x$greedy)) ; x$greedy <-NULL
    if (length(x$mymethod)>0){
      cat('  * User-defined Algorithm Method:\n')
      print(x$mymethod)
    }
    if (x$myterminateflag){
      cat('  * User-defined Termination:\n')
      print(x$terminate)
    }
    if (length(x$output)>0){
      cat('-Customizable output for specialized function:\n')
      print(x$output)
    }
  } else {
    str(x)
  }
}