fix.names<-function(nam, upper=FALSE, lower=FALSE)
  {
   ## fix a set of names to remove blanks
    ##   quote (single and double)
     ##  other bad stuff
  ##   convert all names to upper or lower case

    if(missing(upper)) upper=FALSE
    if(missing(lower)) lower=FALSE
    gnam = nam
    gnam= gsub(" ", "_", gnam)
    gnam=gsub("\'", "_", gnam)
    gnam= gsub("\"", "_", gnam)
    
    if(lower)gnam =  tolower(gnam)
    if(upper) gnam =   toupper(gnam)

    
    
    return(gnam)
    
  }
###  fix.names(sisroster$fullnam)
