# file:    cohensD.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 16 November 2013

# cohensD() computes Cohen's d measure of effect size for the difference 
# between two sample means
cohensD <- function(x = NULL, y = NULL, data = NULL, method = "pooled",  mu = 0, formula=NULL ) {
  
  #### there's only a limited number of meaningful input patterns. ENFORCE this ###
  
  # determine which variables the user has input
  userInput <- !(c( x = missing(x),
                    y = missing(y),
                    data = missing(data),
                    method = missing(method),
                    mu = missing(mu),
                    formula = missing(formula)
                ))
  scenario <- "bad"
  
  # user inputs x 
  if(all(  userInput == c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)  )) { 
    
    if( is(x,"numeric")) { 
      # good case 1: x=numeric -> single sample cohensD with mu=0
      scenario <- "one.sample"
    }
    
    if( is(x,"formula")) {
      # good case 2: x=formula -> two sample cohensD with variables in workspace
      scenario <- "two.sample.formula"
      
      # since no y is specified, create it now
      y <- try( eval( model.frame(formula = x), envir=parent.frame() ), silent=TRUE)
      if( is(y,"try-error") ) { 
        stop( "variables specified in the formula do not exist or are different lengths")
      }
    } 
  }
  
  
  # good case 3: x=numeric, mu=numeric -> single sample cohensD
  if(all(  userInput == c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE)  )) { # user inputs x,mu
    if( is(x,"numeric") ) { # check input type for x
      if( is(mu,"numeric") & length(mu)==1 ) { # check input for mu
        scenario <- "one.sample"
      }
    }
  }
                 
  
  # user inputs x,y
  if(all(  userInput == c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)  )) { 
    
    # good case 4: x=numeric, y=numeric -> two sample cohensD
    if( is(x,"numeric") & is(y,"numeric") ) {
      scenario <- "two.sample.grouped"
    }
    
    # good case 5: x=formula, y=df -> two sample cohensD from formula input
    if( is(x,"formula") & is(y,"data.frame") ) {
      scenario <- "two.sample.formula"
    }
    
  }
    
  # good case 6: x=formula, data=df -> two sample cohensD
  if(all(  userInput == c(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE)  )) { # user inputs x,data
    if( is(x,"formula") & is(data,"data.frame") ) {
      y <- data
      data <- NULL
      scenario <- "two.sample.formula"
    } 
  }  
  
  # good case 7: formula=formula, data=df -> two sample cohensD
  if(all(  userInput == c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)  )) { # user inputs formula,data
    if( is(formula,"formula") & is(data,"data.frame") ) {
      y <- data
      data <- NULL
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
    }  
  }  
  
  # good case 8: x=formula, method=character
  if(all(  userInput == c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE)  )) { 
    
    if( is(x,"formula") & is(method,"character") & length(method)==1 ) {
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
      
      # since no y is specified, create it now
      y <- try( eval( model.frame(formula = x), envir=parent.frame() ), silent=TRUE)
      if( is(y,"try-error") ) { 
        stop( "variables specified in the formula do not exist or are different lengths")
      }
    } 
  }
  
  # user inputs x,y,method
  if(all(  userInput == c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE)  )) { 
    
    # good case 9: x=numeric, y=numeric, method=character -> two sample cohensD
    if( is(x,"numeric") & is(y,"numeric") & is(method,"character") & length(method)==1) {
      scenario <- "two.sample.grouped"
    }
    
    # good case 10: x=formula, y=df, method=character -> two sample cohensD from formula input
    if( is(x,"formula") & is(y,"data.frame") & is(method,"character") & length(method)==1) {
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
    }
    
  }
  
  # good case 11: x=formula, data=df, method=character -> two sample cohensD
  if(all(  userInput == c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE)  )) { # user inputs x,data,method
    if( is(x,"formula") & is(data,"data.frame") & is(method,"character") & length(method)==1 ) {
      y <- data
      data <- NULL
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
    } 
  }    
  
  # good case 12: formula=formula, data=df, method=character -> two sample cohensD
  if(all(  userInput == c(FALSE,FALSE,TRUE,TRUE,FALSE,TRUE)  )) { # user inputs formula,data,method
    if( is(formula,"formula") & is(data,"data.frame") & is(method,"character") & length(method)==1 ) {
      y <- data
      data <- NULL
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
    }  
  }  

  # good case 13: formula=formula, method=character
  if(all(  userInput == c(FALSE,FALSE,FALSE,TRUE,FALSE,TRUE)  )) { 
    
    if( is(formula,"formula") & is(method,"character") & length(method)==1 ) {
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
      
      # since no y is specified, create it now
      y <- try( eval( model.frame(formula = x), envir=parent.frame() ), silent=TRUE)
      if( is(y,"try-error") ) { 
        stop( "variables specified in the formula do not exist or are different lengths")
      }
    } 
  }
  
  # good case 14: formula=formula
  if(all(  userInput == c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)  )) { 
    
    if( is(formula,"formula") ) {
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
      
      # since no y is specified, create it now
      y <- try( eval( model.frame(formula = x), envir=parent.frame() ), silent=TRUE)
      if( is(y,"try-error") ) { 
        stop( "variables specified in the formula do not exist or are different lengths")
      }
    } 
  }
  
  # good case 15: formula=formula, x=df -> two sample cohensD
  if(all(  userInput == c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)  )) { # user inputs formula,x 
    if( is(formula,"formula") & is(x,"data.frame") ) {
      y <- x
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
    } 
  } 
  
  # good case 16: formula=formula, x=df, method=character -> two sample cohensD
  if(all(  userInput == c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE)  )) { # user inputs formula,x,method
    if( is(formula,"formula") & is(x,"data.frame") & is(method,"character") & length(method)==1 ) {
      y <- x
      x <- formula
      formula <- NULL
      scenario <- "two.sample.formula"
      if( method=="paired") { # issue warning about case ordering...
        warning( "calculating paired samples Cohen's d using formula input. Results will be incorrect if cases do not appear in the same order for both levels of the grouping factor") 
      } 
    } 
  } 
  
  ####### throw error if scenario is bad  ####### 
  if( scenario == "bad" ) { 
    stop( "arguments specified do not appear to correspond to a meaningful cohen's d calculation")
  }
  
  
  ####### run one-sample calculation if scenario is one sample  ####### 
  if ( scenario == "one.sample" ) { 
    x <- x[!is.na(x)]
    d <- abs(mean(x) - mu) / sd(x) 
    return( d )
  }  
  
  ####### convert to grouped data if scenario is two sample formula  #######  
  if( scenario == "two.sample.formula" ) {
    
    outcome <- eval( x[[2]], y ) # not super happy about using eval here...
    group <- eval( x[[3]], y )
    group <- as.factor(group) 
    if (nlevels(group) != 2L) {
      stop("grouping factor must have exactly 2 levels")
    }
    x <- split(outcome,group)
    y <- x[[2]]
    x <- x[[1]]
    scenario <- "two.sample.grouped"
    
  }
  
  ######## do a two sample calculation ######## 
  if( scenario == "two.sample.grouped" ) { 
    
    # check method
    if( !( method %in% c("x.sd","y.sd","pooled","corrected","raw","paired","unequal" ))) {
      stop( '"method" must be "x.sd","y.sd","pooled","corrected","raw","paired" or "unequal"')
    }
    
    # remove missing data
    if( method == "paired" ) { # remove complete cases for paired
      if( length(x) != length(y) ) {
        stop( "paired samples cohen's d requires samples of the same size")
      }
      ind <- !is.na(x) & !is.na(y)
      x <- x[ind]      
      y <- y[ind]
    } else { # remove individual cases for non-paired
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]      
    }
    
    
    # function to compute a pooled standard deviation
    pooledSD <- function(x,y,debias = TRUE) {  
      sq.devs <- (c( x-mean(x), y-mean(y) ))^2
      n <- length(sq.devs)
      if(debias) { psd <- sqrt(sum(sq.devs)/(n-2)) }
      else { psd <- sqrt(sum(sq.devs)/n)}
      return(psd)
    }
    
    # calculate cohens d
    mean.diff <- mean(x) - mean(y)
    sd.est <- switch( EXPR = method,             
                      "x.sd" = sd(x),
                      "y.sd" = sd(y),
                      "pooled" = pooledSD(x,y),
                      "corrected" = pooledSD(x,y),
                      "raw" = pooledSD(x,y,FALSE),
                      "paired" = sd(x-y),
                      "unequal" = sqrt( (var(x)+var(y))/2 )
    )               
    d <- mean.diff / sd.est
    if( method == "corrected") { 
      n <- length(x) + length(y)
      d <- d * (n - 3)/(n-2.25)
    }
    
    return(abs(d))
  }
  

  # check
  stop( "how did I get here?")
 
}
