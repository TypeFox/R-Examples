# file:    pairedSamplesTTest.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 30 January 2014


pairedSamplesTTest <- function(
  formula,
  data=NULL,
  id=NULL,
  one.sided = FALSE,
  conf.level=.95
) {
  
  
  
  # check that the user has input a formula
  if( missing(formula) ) { stop( '"formula" argument is missing, with no default')}
  if( !is( formula, "formula")) { stop( '"formula" argument must be a formula')}
  
    
  if( length(formula)==2) { ############ ONE-SIDED FORMULA ############ 
    
    # read off the formula
    vars <- all.vars( formula ) 
    if( length(vars) != 2 ) stop( "one-sided 'formula' must contain exactly two variables" )
    outcome <- vars
    gp.names <- vars
    group <- NA
    id <- NA
    
    # check the data
    if( !missing(data) ) { # is there a data frame?
      
      # it needs to be data frame, because a matrix can't 
      # contain both factors and numeric variables
      if( !is(data,"data.frame") ) stop ( "'data' is not a data frame")
      
      # check that all three variables are in the data frame
      if( !( vars[1] %in% names(data)) ) {
        stop( paste0( "'", vars[1], "' is not the name of a variable in '", deparse(substitute(data)), "'" ))
      }
      if( !( vars[2] %in% names(data)) ) {
        stop( paste0( "'",vars[2],"' is not the name of a variable in '", deparse(substitute(data)), "'" ))
      }      
      
      # truncate the data frame
      data <- data[,vars]
      
    } else {
      
      # check that all variables exist in the workspace
      workspace <- objects( parent.frame())
      
      # check that all three variables are in the data frame
      if( !( vars[1] %in% workspace) ) {
        stop( paste0( "'", vars[1], "' is not the name of a variable in the workspace" ))
      }
      if( !( vars[2] %in% workspace) ) {
        stop( paste0( "'", vars[2],"' is not the name of a variable in the workspace" ))
      }
      
      # copy variables into a data frame if none is specified, and
      # check that the variables are appropriate for a data frame
      ff <- as.formula( paste( "~", vars[1], "+", vars[2]))
      data <- try( eval( model.frame( formula = ff, na.action = na.pass ), 
                         envir=parent.frame() ), silent=TRUE)
      if( is(data,"try-error") ) { 
        stop( "specified variables cannot be coerced to data frame")
      }
      
    }
    
    # check classes of the variables
    if( !is( data[,vars[1]], "numeric" ) ) stop( paste0( "'", vars[1], "' is not numeric" ) ) 
    if( !is( data[,vars[2]], "numeric" ) ) stop( paste0( "'", vars[2], "' is not numeric" ) )  
    
    # remove missing data
    exclude.id <- is.na( data[,vars[1]]) | is.na( data[,vars[2]])
    if( sum(exclude.id) > 0){
      warning( paste( sum(exclude.id), "case(s) removed due to missingness" ) )
    }
    data <- data[!exclude.id,,drop=FALSE]
    
    # create data matrix for later with dummy column 
    WF <- cbind(rep.int(NA,nrow(data)),data)
    
    ############ check the one-sided option ############ 
    if( length(one.sided) !=1 ) stop( "invalid value for 'one.sided'" )
    if( one.sided == FALSE ) { # two sided
      alternative <- "two.sided"
    } else {
      if( one.sided == vars[1] ) { # first variable is the bigger one
        alternative <- "greater"
      } else { 
        if( one.sided == vars[2] ) { # second variable is the bigger one
          alternative <- "less"
        } else {
          stop( "invalid value for 'one.sided'" )
        }
      }
    }
    
    
     
  } else {  ############ TWO-SIDED FORMULA ############  
    
    ############ check formula / id combination ############ 
    
    # check that the user has specified an id that might map onto a variable name
    if( !missing(id) ) { # yes, there's an id...
      
      # is it a character of length one?    
      if( !is(id,"character") | length(id) !=1 ) { 
        stop( '"id" argument does not specify the name of a valid id variable')
      }
      
      # if there's an id, then the formula must be of the form DV ~ IV 
      if( length( formula ) !=3 ) stop( 'invalid value for "formula" argument' )
      vars <- all.vars( formula ) 
      if( length( vars) !=2 ) stop( 'invalid value for "formula" argument' )
      outcome <- vars[1]
      group <- vars[2]
      
      
    } else { # no, there isn't...
      
      # if there's no id, then the formula must specify the id variable in a
      # lme4-like fashion... either this:  DV ~ IV + (id) or DV ~ IV + (1|id).
      # this functionality is not properly, but my own sense of elegance 
      # makes me want to be able to specify the full model via the formula
      
      meets.sneaky.case <- FALSE
      
      if( length( formula)==3 ) { # must be a two sided formula...
        outcome <- all.vars(formula[[2]]) # pull the outcome variable [to be checked later]
        if( length( outcome)==1 ) { # must contain only one outcome variable...
          
          rhs <- formula[[3]] # grab the right hand side of the formula
          if( is( rhs, "call" ) &&  # RHS must be a call
                length( rhs)==3 &&  # must be a binary operation
                deparse( rhs[[1]]) == "+"  # that operation must be +
          ) { 
            terms <- strsplit( deparse(rhs), split="+", fixed=TRUE)[[1]] # split by +
            if( length(terms) == 2 ) { # must have only two terms...
              terms <- gsub(" ","",terms) # deblank
              
              id.candidate <- grep("\\(.*\\)",terms) # id variable must have (.*) in it
              if( length( id.candidate) == 1) { # there can be only 1
                id <- terms[id.candidate] # grab that term
                group <- terms[3-id.candidate] # assume the other one is the group [it is checked later] 
                
                # does it match lme4-like (1|id) ?
                if( length(grep( "^\\(1\\|.*\\)$", id ))==1 ) {
                  id <- gsub( "^\\(1\\|", "", id ) # delete the front bit
                  id <- gsub( "\\)$", "", id ) # delete the back bit
                  formula <- as.formula( paste(outcome, "~", group) ) # truncated formula
                  meets.sneaky.case <- TRUE
                  
                } else {
                  
                  # alternatively, does it match (id) ?
                  if( length(grep( "^\\(.*\\)$", id ))==1  ) {
                    id <- gsub( "^\\(", "", id ) # delete the front bit
                    id <- gsub( "\\)$", "", id ) # delete the back bit
                    formula <- as.formula( paste(outcome, "~", group) ) # truncated formula
                    meets.sneaky.case <- TRUE                 
                    
                  }  
                } 
              }   
            }
          }
        }      
      }
      
      if( !meets.sneaky.case ) stop( "no 'id' variable specified")
      
    }
    
    
    ############ check data frame ############ 
    
    # at this point we know that outcome, vars and id are all character
    # vectors that are supposed to map onto variables either in the workspace
    # or the data frame
    
    # if the user has specified 'data', check that it is a data frame that 
    # contains the outcome, group and id variables.
    if( !missing(data) ) { 
      
      # it needs to be data frame, because a matrix can't 
      # contain both factors and numeric variables
      if( !is(data,"data.frame") ) stop ( "'data' is not a data frame")
      
      # check that all three variables are in the data frame
      if( !( outcome %in% names(data)) ) {
        stop( paste0( "'", outcome, "' is not the name of a variable in '", deparse(substitute(data)), "'" ))
      }
      if( !( group %in% names(data)) ) {
        stop( paste0( "'",group,"' is not the name of a variable in '", deparse(substitute(data)), "'" ))
      }
      if( !( id %in% names(data)) ) {
        stop( paste0( "'",id,"' is not the name of a variable in '", deparse(substitute(data)), "'" ))
      }
      
      
    } else {
      
      # check that all variables exist in the workspace
      workspace <- objects( parent.frame())
      
      # check that all three variables are in the data frame
      if( !( outcome %in% workspace) ) {
        stop( paste0( "'", outcome, "' is not the name of a variable in the workspace" ))
      }
      if( !( group %in% workspace) ) {
        stop( paste0( "'",group,"' is not the name of a variable in the workspace" ))
      }
      if( !( id %in% workspace) ) {
        stop( paste0( "'",id,"' is not the name of a variable in the workspace" ))
      }
      
      # copy variables into a data frame if none is specified, and
      # check that the variables are appropriate for a data frame
      ff <- as.formula( paste( outcome, "~", group, "+", id))
      data <- try( eval( model.frame( formula = ff, na.action = na.pass ), 
                         envir=parent.frame() ), silent=TRUE)
      if( is(data,"try-error") ) { 
        stop( "specified variables cannot be coerced to data frame")
      }
      
    }
    
    # subset the data frame
    data <- data[, c(outcome,group,id) ]
    
    
    ############ check classes for outcome, group and id ############ 
    
    # at this point we have a data frame that is known to contain 
    # outcome, group and id. Now check that they are of the appropriate
    # type to run a t-test
    
    # outcome must be numeric  
    if( !is(data[,outcome],"numeric") ) stop( "outcome variable must be numeric")
    
    # group should be a factor with two-levels. issue warnings if it only
    # has two unique values but isn't a factor, or is a factor with more 
    # than two levels but only uses two of them. 
    
    if( is(data[,group], "factor") ) { # it's a factor
      
      if( nlevels( data[,group]) <2 ) { # fewer than two levels
        stop( "grouping variable does not have two distinct levels")
      } 
      
      if( nlevels( data[,group]) >2 ) { # more than two levels
        if( length( unique( data[,group] ))==2 ) { # but only two of them are used...
          warning( "grouping variable has unused factor levels")
          data[,group] <- droplevels( data[,group]) 
          
        } else { # too many levels in use
          stop( "grouping variable has more than two distinct values")
        }
      }
      
    } else { # it's not a factor
      
      if( length( unique( data[,group] ))==2 ) { # if it happens to have 2 unique values...
        warning( "group variable is not a factor" ) # warn the user
        data[,group] <- as.factor( data[,group]) # coerce and continue...
        
      } else {
        stop( "grouping variable must contain only two unique values (and should be a factor)")
      }
      
    }
    
    
    # id should be a factor. issue a warning if it isn't
    if( !is( data[,id], "factor" )) warning( "id variable is not a factor")
    
    
    
    ############ check the one-sided option ############ 
    
    # group names
    gp.names <- levels(data[,group])
    
    # check alternative
    if( length(one.sided) !=1 ) stop( "invalid value for 'one.sided'" )
    if( one.sided == FALSE ) { # two sided
      alternative <- "two.sided"
    } else {
      if( one.sided == gp.names[1] ) { # first factor level
        alternative <- "greater"
      } else { 
        if( one.sided == gp.names[2] ) { # second factor level
          alternative <- "less"
        } else {
          stop( "invalid value for 'one.sided'" )
        }
      }
    }
    
    
    ############ check cases and restructure ############ 
    
    # check that we have the right number of cases?
    tt <- table( data[,id], data[,group] )
    if( any(tt > 1) ) stop( "there too many observations for some cases" )
    
    # find cases to remove
    exclude.id <- tt[,1] !=1 | tt[,2] != 1   # exclude if the relevant row is missing...
    more.bad <- as.character(unique(data[apply( is.na(data[,c(outcome,group)]), 1, any ),id])) # or if it has NA
    exclude.id[ more.bad ] <- TRUE   
    if( sum(exclude.id) > 0){
      warning( paste( sum(exclude.id), "case(s) removed due to missingness" ) )
    }
    exclude.id <- rownames(tt)[exclude.id]
    
    # remove bad cases if they exist
    bad.cases <- data[,id] %in% exclude.id
    data <- data[ !bad.cases,, drop=FALSE ]
    
    # Convert to wide form for later
    WF <- longToWide( data, formula )
    
  }
  
  
  ############ check the confidence level ############ 
  
  if( !is(conf.level,"numeric") |
        length( conf.level) != 1 |
        conf.level < 0 |
        conf.level > 1  
  ) {
    stop( '"conf.level" must be a number between 0 and 1')    
  }
  
  ############ do the statistical calculations ############ 
  
  # pass to t.test
  htest <- t.test( x = WF[,2], y=WF[,3], paired=TRUE,
                   alternative=alternative, conf.level=conf.level )
  
  # cohens D
  d <- cohensD( x= WF[,2], y=WF[,3], method="paired" )
  
  # descriptives
  gp.means <- sapply( WF[,-1], mean )
  gp.sd <- sapply( WF[,-1], sd )
  
  # add the difference scores to the descriptives
  gp.means <- c( gp.means, mean(WF[,2]-WF[,3]))
  gp.sd <- c( gp.sd, sd(WF[,3]-WF[,2]))
  
  ############ output ############ 
  
  # create output structure
  TT <- list( 
    t.statistic = htest$statistic,
    df = htest$parameter,
    p.value = htest$p.value,
    conf.int = htest$conf.int,
    conf = conf.level,
    mean = gp.means, 
    sd = gp.sd, 
    outcome = outcome,
    group = group,
    group.names = gp.names,
    id = id,
    mu = NULL,
    alternative = alternative,
    method = "Paired samples t-test",
    effect.size = d
  )
  
  # specify the class and return
  class(TT) <- "TTest"
  return(TT)
  
  
}


