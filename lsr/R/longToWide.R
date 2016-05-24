# file:    longToWide.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 29 June 2013

# longToWide() takes a long-form data frame and converts it to a wide-form data frame. 
# Like it's companion function wideToLong() it's not as flexible as cast() and melt(),
# but it is easier to use.
longToWide <- function( data, formula, sep="_") {
  
  within <- all.vars(formula[-2])
  v.names <- all.vars(formula[-3])
  idvar <- setdiff(names(data),c(within,v.names)) 
  
  if( length(within)>1 ) { 
    collapsed.treatments <- apply(as.matrix(data[,within]),1,paste,collapse=sep)
    data <- data[,setdiff(names(data),within)] # delete split treatments
    data$within <- collapsed.treatments # append collapsed treatment
    within <- "within"
  }
  times <- unique( data[,within]) # measure 'time' names
  varying <- list()
  for( i in seq_along(v.names) ) varying[[i]] <- paste(v.names[i],times, sep=sep)
  
  x<-reshape( data, idvar=idvar, varying=varying, direction="wide", times=times, v.names=v.names, timevar=within)
  rownames(x) <- NULL
  return(x)
  
}