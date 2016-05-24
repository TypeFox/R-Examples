# file:    associationTest.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 23 January 2014



associationTest <- function( formula, data=NULL ) {
  
  ############ check  formula ############ 
  
  # check that the user has input a formula
  if( missing(formula) ) { stop( '"formula" argument is missing, with no default')}
  if( !is( formula, "formula")) { stop( '"formula" argument must be a formula')}
  
  # the formula must be of the form ~ VAR1 + VAR2
  if( length( formula ) !=2 ) stop( '"formula" argument must be one-sided with two variables' )
  vars <- all.vars( formula ) 
  if( length( vars) !=2 ) stop( '"formula" argument must be one-sided with two variables' )
  
  ############ check data frame ############ 
  
  if( !missing(data) ) { 
    
    # it needs to be data frame, because a matrix can't 
    # contain actors
    if( !is(data,"data.frame") ) stop ( "'data' argument must be a data frame")
    
    # check that both variables are in the data frame
    if( !( vars[1] %in% names(data)) ) {
      stop( paste0( "'", vars[1], "' is not the name of a variable in '", deparse(substitute(data)), "'" ))
    }
    if( !( vars[2] %in% names(data))) {
      stop( paste0( "'", vars[2], "' is not the name of a variable in '", deparse(substitute(data)), "'" ))
    }
    
  } else {
    
    # check that all variables exist in the workspace
    workspace <- objects( parent.frame())
    
    # check that both variables are in the workspace
    if( !( vars[1] %in% workspace)) {
      stop( paste0( "'", vars[1], "' is not the name of a variable in the workspace" ))
    }
    if( !( vars[2] %in% workspace)) {
      stop( paste0( "'", vars[2], "' is not the name of a variable in the workspace" ))
    }
    
    # copy variables into a data frame if none is specified, and
    # check that the variables are appropriate for a data frame.
    # need to retain the missing values for later
    data <- try( eval( model.frame( formula = formula, na.action = na.pass ), 
                       envir=parent.frame() ), silent=TRUE )
    if( is(data,"try-error") ) { 
      stop( "specified variables cannot be coerced to data frame")
    }
  }
  
  # subset the data frame
  data <- data[, vars ]
  
  ############ check data frame ############   
  
  # check that both variables are factors
  if( !is(data[,vars[1]],"factor")) {
    stop( paste0( "'",vars[1],"' is not a factor"))
  }
  if( !is(data[,vars[2]],"factor")) {
    stop( paste0( "'",vars[2],"' is not a factor"))
  }    
    

  
  # check for missing data & print warning if needed
  missingData <- is.na( data[,vars[1]]) | is.na( data[,vars[2]])
  data <- data[!missingData,] 
  if(sum(missingData) > 0) {
    warning( paste(sum(missingData)), " case(s) removed due to missingness" )
  }
  
  
  ############ do the test ############   
  
  # get the variable names  
  f <- xtabs( formula, data )
    
  # run the corresponding chi-square test, suppressing the warning,
  # replacing it with our own if it exists
  old.warn <- options(warn=2) # convert warnings to errors
  htest <- try( chisq.test( x=f ), silent=TRUE ) # try the test
  need.warning <- FALSE # assume warning
  
  # check for failures:
  if( class(htest) == "try-error" ) {
    
    # handle the case when the "error" was the warning we're catching
    need.warning <- length( grep("Chi-squared approximation may be incorrect",htest )) > 0
    if( need.warning ) { # if it was a ties problem... 
      options(warn=-1) # suppress warnings completely for the next run...
    } else { # if not...
      options( old.warn ) # reset warning state and let R throw what it likes...
    }
    
    # now run the test again & compute effect size
    htest <- chisq.test( x=f )
  } 
  effectsize <- cramersV( f )
  options( old.warn )   # reset warnings 
  
  # make a note of whether yates' correction has been applied
  yates <- FALSE
  if( nrow(f)==2 & ncol(f)==2 ) yates <- TRUE
  
  
  ############ output ############ 
  
  # create output structure
  out <- list( 
    variables = vars,
    observed = htest$observed,
    expected = htest$expected,
    difference = htest$observed - htest$expected,
    statistic = htest$statistic,
    df = htest$parameter,
    p.value = htest$p.value,
    effect.size = effectsize,
    yates = yates,
    warn = need.warning
  )
  class( out ) <- "assocTest"
  
  # throw the warning
  if( need.warning ) { 
    warning( "Expected frequencies too small: chi-squared approximation may be incorrect")
  }
  
  return(out)
  
}


# print method
print.assocTest <- function(x, ...) {
  
  nDigits <- 3
  
  # print the name of the test
  cat("\n     Chi-square test of categorical association\n\n")
  
  # print the data variable
  cat( "Variables:  ", paste(x$variables,collapse=", "), "\n\n" )
  
  # print the hypotheses being tested
  cat( "Hypotheses: \n")
  cat( "   null:        variables are independent of one another\n")
  cat( "   alternative: some contingency exists between variables\n")
  cat("\n")
  
  # print out the observed
  cat( "Observed contingency table:\n")
  print( x$observed)
  cat("\n")
  
  # print out the expected
  cat( "Expected contingency table under the null hypothesis:\n")
  print( x$expected, digits=nDigits)
  cat("\n")
  
  # print test statistics
  cat( "Test results: \n")
  cat( "   X-squared statistic: ", round( x$statistic, nDigits ), "\n" )
  cat( "   degrees of freedom: ", round( x$df, nDigits ), "\n" )
  pp <- ifelse( x$p.value < .001, "<.001", round( x$p.value, nDigits ))
  cat( "   p-value: ", pp, "\n")
  cat( "\n")
  
  # print other things
  cat( "Other information: \n") 
  cat( "   estimated effect size (Cramer's v): ", round( x$effect.size, nDigits), "\n" )  
  if( x$yates) cat( "   Yates' continuity correction has been applied\n")
  if( x$warn) cat( "   warning: expected frequencies too small, results may be inaccurate\n")
  cat( "\n")

  
  return( invisible(x) )
  
} 

