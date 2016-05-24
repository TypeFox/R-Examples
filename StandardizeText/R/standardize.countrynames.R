standardize.countrynames <-
function(input,input.column=NULL,standard="default",standard.column=NULL,only.names=FALSE,na.rm=FALSE,suggest="prompt",print.changes=TRUE,verbose=FALSE) {
  #Initialize country data
  temp.env <- new.env()
  data(country.names,envir=temp.env)
  data(country.regex,envir=temp.env)
  country.names <- get("country.names",envir=temp.env)
  country.regex <- get("country.regex",envir=temp.env)
  #Check for use of predefined country names
  if(check.value(standard,"string")) {
    standard.names <- switch(standard, default = country.names[[1]],
            imf = country.names[[2]], iso = country.names[[3]], 
            pwt = country.names[[4]], wb = country.names[[5]], 
            who = country.names[[6]], NULL)
    if (is.null(standard.names)) {
      cat(sprintf("Error in 'standard': '%s' not a recognized name set\n",standard))
      return()
    } else {
      standard <- standard.names
    }
  }
  return(standardize.names(input,input.column,standard,standard.column,regex=country.regex[,1],codes=country.regex[,2],match=FALSE,only.names,na.rm,suggest,print.changes,verbose))
}
