# file:    printTTest.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 30 January 2014


print.TTest <- function( x, ... ) {
  
  # function to force equal digit printing
  makeTxt <- function(x, nDigits, naPrint ) {
    n <- dim(x)
    format <- paste0("%.",nDigits,"f")
    txt <- sprintf( format, x)
    txt <- gsub("NA", naPrint, txt, fixed=TRUE)
    txt <- matrix(txt,n[1],n[2], dimnames=dimnames(x))
    return(txt)
  }
  
  # function to print the text matrix
  printTxt <- function( txt ) {
    print.default( txt, quote=FALSE, right=TRUE)
  }
  
  # number of digits to round to
  nDigits <- 3
  
  # compute means and standard deviations
  descriptives <- rbind( x$mean, x$sd )
  rownames(descriptives) <- c("   mean","   std dev.")
  
  # print test type
  cat( "\n  ", x$method, "\n\n" )
  
  # which variables are being tested
  if( x$method == "One sample t-test") {
    cat( "Data variable:  ", x$outcome, "\n" )
    colnames( descriptives ) <- gsub( "^.*\\$","", x$outcome )
  } 
  if( x$method == "Student's independent samples t-test" | 
      x$method == "Welch's independent samples t-test") {
    cat( "Outcome variable:  ", x$outcome, "\n" )
    cat( "Grouping variable: ", x$group, "\n" )         
  }
  if( x$method == "Paired samples t-test" ) {
    if( !is.na( x$id) ) { # two sided formula...
      cat( "Outcome variable:  ", x$outcome, "\n" )
      cat( "Grouping variable: ", x$group, "\n" )   
      cat( "ID variable:       ", x$id, "\n" )
    } else {
      cat( "Variables: ", x$outcome[1], ",", x$outcome[2], "\n" )
    }
    colnames( descriptives ) <- c(x$group.names, "difference")
  }
  
  cat("\n")
  
  # print the descriptives
  cat( "Descriptive statistics: \n")
  descriptives.txt <- makeTxt( descriptives, nDigits, "NA" ) # textify
  printTxt( descriptives.txt )
  cat("\n")
  
  # print the hypotheses being tested
  cat( "Hypotheses: \n")
  if( x$method=="One sample t-test") { # one sample null...
    
    # two-sided test
    if( x$alternative == "two.sided" ) {
      cat( "   null:        population mean equals", x$mu, "\n" )
      cat( "   alternative: population mean not equal to", x$mu, "\n" )
    }
    
    # greater-than test
    if( x$alternative == "greater" ) {
      cat( "   null:        population mean less than or equal to", x$mu, "\n")
      cat( "   alternative: population mean greater than", x$mu, "\n" )
    }
    
    # less-than test
    if( x$alternative == "less" ) {
      cat( "   null:        population mean greater than or equal to", x$mu, "\n" )
      cat( "   alternative: population mean less than", x$mu, "\n" )
    }
    
  } else {
     if( x$method=="Paired samples t-test" ) { # paired sample null...
     
      # two-sided test
      if( x$alternative == "two.sided" ) {
        cat( "   null:        population means equal for both measurements\n" )
        cat( "   alternative: different population means for each measurement\n" )
      }
      
      # greater than test
      if( x$alternative == "greater" ) {
        cat( "   null:        population means are equal, or smaller for measurement",paste0("'",x$group.names[1],"'"),"\n" )
        cat( "   alternative: population mean is larger for measurement",paste0("'",x$group.names[1],"'"),"\n" )
      }
      
      # less than test
      if( x$alternative == "less" ) {
        cat( "   null:        population means are equal, or smaller for measurement", paste0("'",x$group.names[2],"'"),"\n" )
        cat( "   alternative: population mean is larger for measurement",paste0("'",x$group.names[2],"'"),"\n" )
      }
      
      
    } else { # two samples null...
      
      # two-sided test
      if( x$alternative == "two.sided" ) {
        cat( "   null:        population means equal for both groups\n" )
        cat( "   alternative: different population means in each group\n" )
      }
      
      # greater than test
      if( x$alternative == "greater" ) {
        cat( "   null:        population means are equal, or smaller for group",paste0("'",x$group.names[1],"'"),"\n" )
        cat( "   alternative: population mean is larger for group",paste0("'",x$group.names[1],"'"),"\n" )
      }
      
      # less than test
      if( x$alternative == "less" ) {
        cat( "   null:        population means are equal, or smaller for group", paste0("'",x$group.names[2],"'"),"\n" )
        cat( "   alternative: population mean is larger for group",paste0("'",x$group.names[2],"'"),"\n" )
      }
      
      
    }
  }
  cat("\n")
  
  # inferential statistics
  cat( "Test results: \n")
  cat( "   t-statistic: ", round( x$t.statistic, nDigits ), "\n" )
  cat( "   degrees of freedom: ", round( x$df, nDigits ), "\n" )
  pp <- ifelse( x$p.value < .001, "<.001", round( x$p.value, nDigits ))
  cat( "   p-value: ", pp, "\n")
  cat( "\n")
  
  # inferential statistics
  cat( "Other information: \n") 
  ci.str <- round( x$conf.int, nDigits )
  ci.str <- paste0( round(x$conf*100), "% confidence interval:  [", ci.str[1], ", ", ci.str[2], "]" )
  if( x$alternative == "two.sided" ) {
    ci.str <- paste( "two-sided", ci.str )
  } else {
    ci.str <- paste( "one-sided", ci.str )
  }
  cat( "  ",ci.str, "\n" )
  cat( "   estimated effect size (Cohen's d): ", round( x$effect.size, nDigits), "\n" )  
  cat( "\n")
  
}