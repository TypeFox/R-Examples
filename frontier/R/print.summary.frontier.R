print.summary.frontier <- function( x, effic = x$printEffic, ... ) {

   if( x$modelType == 1 ) {
      cat( "Error Components Frontier (see Battese & Coelli 1992)\n" )
   } else if( x$modelType == 2 ) {
      cat( "Efficiency Effects Frontier (see Battese & Coelli 1995)\n" )
   } else {
      stop( "unknown model type '", x$modelType, "'" )
   }
   if( x$ineffDecrease ) {
      cat( "Inefficiency decreases the endogenous variable",
	"(as in a production function)\n" )
   } else {
      cat( "Inefficiency increases the endogenous variable",
	"(as in a cost function)\n" )
   }
   if( x$logDepVar == 1 ) {
      cat( "The dependent variable is logged\n" )
   } else {
      cat( "The dependent variable is not logged\n" )
   }
   if( x$nIter == 1 ) {
      if( x$maxit == 1 ) {
         cat( "did not iterate as argument 'maxiter' was set to '1'\n" )
      } else if( x$code == 5 ) {
         cat( "iteration failed:\n",
            "cannot find a parameter vector that results in a log-likelihood value\n",
            "larger than the log-likelihood value of the initial parameters\n",
            sep = "" )
      } else {
         cat( "iteration failed\n" )
      }
   } else {
      cat( "Iterative ML estimation terminated after", x$nIter, "iterations:\n" )
      if( x$code == 1 ) {
         cat( "log likelihood values and parameters of two successive iterations\n",
            "are within the tolerance limit\n", sep = "" )
      } else if( x$code == 5 ) {
         cat( "cannot find a parameter vector that results in a log-likelihood value\n",
            "larger than the log-likelihood value obtained in the previous step\n",
            sep = "" )
      } else if( x$code == 6 ) {
         cat( "search failed on gradient step\n" )
      } else if( x$code == 10 ) {
         cat( "maximum number of iterations reached\n" )
      } else {
         cat( "unknown return code:", x$code, "\n" )
      }
   }
   if( x$nRestart != 0 ) {
      cat( "Multiplied the initial values", x$nRestart, "time(s) by",
         x$restartFactor, "before the search procedure could start\n" )
      cat( "You could try to use different starting values or",
         "try to reduce the step size specified in argument 'searchStep'\n" )
   }

   cat( "\nfinal maximum likelihood estimates\n" )
   printCoefmat( coef( x ), ... )
   cat( "log likelihood value:", x$mleLogl, "\n" )

   if( x$nt == 1 ) {
      cat( "\ncross-sectional data\n" )
      cat( "total number of observations =", x$nob, "\n" )
   } else {
      cat( "\npanel data\n" )
      cat( "number of cross-sections =", x$nn, "\n" )
      cat( "number of time periods =", x$nt, "\n" )
      cat( "total number of observations =", x$nob, "\n" )
      cat( "thus there are", x$nn * x$nt - x$nob,
         "observations not in the panel\n" )
   }

   if( effic ){
      cat( "\nefficiency estimates\n" )
      print( x$effic )
   }

   if( !is.null( x$efficYearMeans ) ) {
      cat( "\nmean efficiency of each year\n" )
      print( x$efficYearMeans )
   }

   cat( "\nmean efficiency:", x$efficMean, "\n" )

   invisible( x )
}