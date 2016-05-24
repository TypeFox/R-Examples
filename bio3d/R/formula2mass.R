"formula2mass" <-
  function(form, sum.mass=TRUE) {
    
    errmsg <- paste("error while parsing formula. \n",
                    " provide input on the form: 'C3 H5 N O1'")
    
    if(class(form)!="character" || missing(form))
      stop(errmsg)
    
    eles <- unlist(strsplit(form, " "))
    
    mass <- c()
    for ( i in 1:length(eles) ) {
      ele <- unlist(strsplit(eles[i], ""))
      
      inds.a <- grep("[0-9]", ele)
      inds.b <- grep("[A-z]", ele)
      
      num <- paste(ele[inds.a], collapse="")
      char <- paste(ele[inds.b], collapse="")
    
      if(nchar(num)==0)
        num <- 1
      if(length(inds.a)==0)
        inds.a <- 0
      
      if(length(inds.b)==0)
        stop(errmsg)
      if(nrow(bounds(inds.a))>1 || nrow(bounds(inds.b))>1 )
        stop(errmsg)
      
      mass <- c(mass, atom2mass(char) * as.numeric(num))
    }
    
    if(sum.mass)
      mass <- sum(mass)
    
    return(mass)
  }
