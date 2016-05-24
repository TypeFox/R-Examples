# file:    importList.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# importList() copies each element of a list into a separate variable in the workspace. 
importList <- function(x, ask = TRUE ) {
  
  if( !is(x,"list") & !is(x,"data.frame")) stop( '"x" must be a list or data frame')
  if( !is(ask,"logical") | length(ask) !=1 ) {
    stop( '"ask" must be a single logical value')
  }
  
  envir = parent.frame() # import to parent environment
  
  vars <- names(x)  # get variable names
  vars <- make.unique(vars) # make sure the names are unique
  vars <- make.names(vars) # convert to legitimate R names 
  
  if( ask ) {
    cat("Names of variables to be created:\n")
    print(vars)
    ans <- NA
    while( ! (ans %in% c("y","n") ) ) {
      ans <- readline("Create these variables? [y/n] ")
    }
    if (ans == "n") { 
      return( invisible(0) )
    }  
  }
  
  for (v in seq_along(vars)) { # for each variable in x:
    assign(x = vars[v], value = x[[v]], envir = envir)
  }
  return( invisible(1) )
  
} 