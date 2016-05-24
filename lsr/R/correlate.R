# file:    correlate.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 26 August 2014

# compute correlation matrix:
#  - automatically removes non-numeric variables from the data
#  - automatically does "pairwise.complete.obs" for handling missing data
#  - can report the result of correlation tests if requested
#  - can do pearson, spearman and kendall
#  - defaults to no test, and to pearson correlation

correlate <- function( x, y=NULL, test=FALSE, corr.method="pearson", p.adjust.method="holm" ) {
  
  # did the user specify two input matrices?
  two.inputs <- !is.null(y)
  
  # allow numeric vectors
  if( is.vector(x) & is.numeric(x) ) {
    x <- data.frame(x)
    call <- match.call()
    n <- call[[2]]
    names(x) <- "x.var"
    if(class(n) == "name") names(x) <- as.character(n)
  }
  if( two.inputs & is.vector(y) & is.numeric(y) ) {
    y <- data.frame(y)
    call <- match.call()
    n <- call[[3]]
    names(y) <- "y.var"
    if(class(n) == "name") names(y) <- as.character(n)
  }
  
  # check input matrices are data frames
  if( !(class(x) %in% c("matrix","data.frame") )) {
    stop( 'x must be a matrix or data frame' ) 
  } 
  if( two.inputs & !(class(y) %in% c("matrix","data.frame") )) {
    stop( 'y must be a matrix or data frame' ) 
  } 
  
  # warn if x and y appear to have the same variables!
  
  # coerce to data frame
  if( class(x) =="matrix" ) x <- as.data.frame(x)
  if( two.inputs & class(y) == "matrix") y <- as.data.frame(y)
  
  # store other args
  args <- c( two.inputs=two.inputs, test=test, 
             corr.method=corr.method, p.adjust.method=p.adjust.method )
  
  # define a function to run correlation test and trap warning message
  getCT <- function( x,y, method ) {
    
    # get the correlation test, trapping the ties problem warning as needed...
    old.warn <- options(warn=2) # convert warnings to errors
    ct <- try( cor.test( x, y, method=method), silent=TRUE ) # try the correlation
    tp <- FALSE # assume no ties problem unless...
    
    # check for failures:
    if( class(ct) == "try-error" ) {
      
      # handle the case when the "error" was a ties warning
      tp <- length( grep("exact p-value with ties",ct )) > 0
      if( tp ) { # if it was a ties problem... 
        options(warn=-1) # suppress warnings completely for the next run...
      } else { # if not...
        options( old.warn ) # reset warning state and let R throw what it likes...
      }
      
      # now run the test again...
      ct <- cor.test( x, y, method=method )
    } 
    options( old.warn ) # reset warnings to original state
    
    return( list( ct=ct, tp=tp) )
    
  }
  
  
  if( !two.inputs ) {
    #### one matrix case ###
  
    # drop categorical variables 
    classes <- sapply(x,"class") # get the classes
    inds <- which( classes %in% c("integer","numeric")) # retain only int and num 
    n.vars <- dim(x)[2] # number of variables
    n.cont <- length(inds) # number of continuous variables
    
    # initialise the output list with empty matrices
    R <- list( correlation = matrix( NA, n.vars, n.vars) )
    rownames( R$correlation ) <- colnames( R$correlation ) <- colnames(x)
    R$sample.size <- R$p.value <- R$correlation  
    R$args <- args
    R$tiesProblem <- FALSE
    
    # run pairwise tests (inefficient looping!)
    for( a in 1:(n.cont-1) ){
      for( b in (a+1):n.cont) {     
        i <- inds[a] # the a-th continuous variable
        j <- inds[b] # the b-th continuous variable
             
       # ct <- cor.test( x[,i], x[,j], method=corr.method ) 
        cttp <- getCT( x[,i], x[,j], corr.method )
        
        # store the output
        R$tiesProblem <- R$tiesProblem | cttp$tp
        R$correlation[j,i] <- R$correlation[i,j] <- cttp$ct$estimate
        R$p.value[j,i] <- R$p.value[i,j] <- cttp$ct$p.value
      }
    }
    
    # adjust p (inefficient duplication here)
    upper.inds <- upper.tri(R$p.value) & !is.na(R$p.value) # cells containing p-values, upper
    R$p.value[upper.inds] <- p.adjust( R$p.value[upper.inds], method=p.adjust.method )
    lower.inds <- lower.tri(R$p.value) & !is.na(R$p.value) # cells containing p-values, lower
    R$p.value[lower.inds] <- p.adjust( R$p.value[lower.inds], method=p.adjust.method )
    
    
    # fill in sample sizes for individual variables
    for( i in 1:n.vars ) {
      R$sample.size[i,i] <- sum(!is.na(x[,i]))
    } 
    
    # and pairwise sample sizes for all variables
    for( i in 1:(n.vars-1) ) {
      for( j in (i+1):n.vars) {
        R$sample.size[j,i] <- R$sample.size[i,j] <- sum(!(is.na(x[,i]) | is.na(x[,j])))
      }
    }
    
  } else {
    #### two matrix case ###   
    
    # track only the categorical variables 
    inds.x <- which(sapply(x,"class") %in% c("integer","numeric")) 
    inds.y <- which(sapply(y,"class") %in% c("integer","numeric"))     
    
    # variable numbers
    n.vars.x <- dim(x)[2] # number of variables in x
    n.vars.y <- dim(y)[2] # number of variables in y
    n.cont.x <- length(inds.x) # number of continuous variables in x
    n.cont.y <- length(inds.y) # number of continuous variables in y
    
    # initialise the output list
    R <- list( correlation = matrix( NA, n.vars.x, n.vars.y) )
    rownames( R$correlation ) <- colnames(x)
    colnames( R$correlation ) <- colnames(y)
    R$sample.size <- R$p.value <- R$correlation  
    R$args <- args
    R$tiesProblem <- FALSE
    
    # run pairwise tests (inefficient looping!)
    for( a in 1:n.cont.x ){
      for( b in 1:n.cont.y) {     
        i <- inds.x[a] # the a-th continuous variable in x
        j <- inds.y[b] # the b-th continuous variable in y
        
        #ct <- cor.test( x[,i], y[,j], method=corr.method ) # run the test
        cttp <- getCT( x[,i], y[,j], corr.method )
        
        
        # store the output
        R$tiesProblem <- R$tiesProblem | cttp$tp
        R$correlation[i,j] <- cttp$ct$estimate
        R$p.value[i,j] <- cttp$ct$p.value
        
        # store sample size
        R$sample.size[i,j] <- sum(!(is.na(x[,i]) | is.na(y[,j]))) 
      }
    }
    
    # adjust p
    R$p.value <- matrix( p.adjust( R$p.value, method=p.adjust.method ), 
                         n.vars.x, n.vars.y,
                         dimnames=dimnames(R$p.value) ) 
    
  }
  
  # define an S3 class
  class(R) <- c("correlate","list")
  return(R)

}


# print method
print.correlate <- function( x, ... ){ 
  
  # fixed properties that should probably be converted to
  # input arguments?  
  nDigits <- 3
  naPrint <- "."
  
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
  
  # print the correlations
  cat("\n")
  cat("CORRELATIONS\n")
  cat("============\n")
  cat("- correlation type: ", x$arg["corr.method"], "\n")
  cat("- correlations shown only when both variables are numeric\n")
  cat("\n")
  
  if( options()$show.signif.stars & x$arg["test"]==TRUE ) { # if significance stars needed...
    
    # function to return the significance stars string
    getSigString <- function(p) {
      if( is.na(p) | p > .1 ) return( "   " )
      if( p > .05 ) return( ".  " )
      if( p > .01 ) return( "*  " ) 
      if( p > .001 ) return( "**" )
      return( "***" )
    }
    
    # function to generate interleaved indices
    interleave <- function(n){
      ord <- vector()
      ord[ seq(1,2*n-1,2) ] <- 1:n
      ord[ seq(2,2*n,2)] <- (n+1):(2*n)
      return(ord)
    }
    
    # print correlation matrix with sig stars
    n <- dim( x$correlation )[2] # number of columns
    txt <- makeTxt( x$correlation, nDigits, naPrint ) # text form of the correlation matrix
    for( i in 1:n ) txt[,i] <- paste0(txt[,i], sapply(x$p.value[,i], getSigString)) # append stars
    colnames(txt) <- paste0( colnames(x$correlation), rep.int("   ",n)) # column names
    rownames(txt) <- rownames( x$correlation ) # row names
    printTxt(txt) 
    
    # print the key
    cat("\n---\n")
    cat("Signif. codes: . = p < .1, * = p<.05, ** = p<.01, *** = p<.001\n")
    
  } else { # if no significance stars needed...
  
    txt <- makeTxt( x$correlation, nDigits, naPrint )
    printTxt( txt )
  }
  
  if( x$arg["test"]==TRUE ) {
  
    # print the p.values
    cat("\n\np-VALUES\n")
    cat("========\n")
    if(  x$args["two.inputs"] ) { nTests <- sum( !is.na( x$correlation))
    } else { nTests <- sum( !is.na( x$correlation)) / 2 }
    cat("- total number of tests run: ", nTests, "\n")
    cat("- correction for multiple testing: ", x$arg["p.adjust.method"],"\n")
    if( x$tiesProblem ) {
      cat("- WARNING: cannot compute exact p-values with ties\n" )
    }
    cat("\n")
    
    # print p-values, forcing nDigits 
    txt <- makeTxt( x$p.value, nDigits, naPrint )
    printTxt( txt )

    
    # print the sample sizes
    cat("\n\nSAMPLE SIZES\n")
    cat("============\n")
    cat("\n")
    txt <- makeTxt( x$sample.size, nDigits=0, naPrint )
    printTxt( txt )

  
  }
  
  # invisbly return the input
  return( invisible(x) )
}



