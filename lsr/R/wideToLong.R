# file:    wideToLong.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 26 September 2013

# wideToLong() takes a wide-form data frame and converts it to a long-form data frame. 
# Like it's companion function longToWide() it's not as flexible as cast() and melt(),
# because it relies on variable names following a very specific naming scheme, but as
# long as the user does adhere to this scheme it's a much easier function to use.
wideToLong <- function( data, within="within", sep="_", split=TRUE) {
  
  ind <- grep(sep,names(data),fixed=TRUE) # indices of variables that are repeated
  idvar <- names(data)[-ind] # names of id variables
  
  # make sure that the id variables do uniquely specify cases
  n.profiles <- dim( unique( as.matrix( data[,idvar,drop=FALSE] ), margin=1 ) )[1] # number of unique id-var profiles
  if( n.profiles < dim(data)[1] ) { # if id variables don't uniquely specify cases
    warning( "Between-subject variables must uniquely specify cases: a case 'id' variable has been added")
    id <- 1:dim(data)[1]
    data <- cbind(data,id)
    names(data) <- make.unique(names(data))
    ind <- grep(sep,names(data),fixed=TRUE) # indices of variables that are repeated
    idvar <- names(data)[-ind] # names of id variables
  }
  
  tmp <- t(as.data.frame(strsplit( names(data[ind]), sep, fixed=TRUE ))) # matrix with split var names
  v.names <- unique(tmp[,1]) # grab the measure var names
  times <- unique(apply( tmp[,-1,drop=FALSE],1,paste,collapse=sep)) # measure 'time' names
  varying <- list()
  for( i in seq_along(v.names) ) varying[[i]] <- names(data)[ind][tmp[,1]==v.names[i]] 
  
  tmp <-make.unique(c(names(data),"withintmp"))
  within.tmp <- tmp[length(tmp)]
  
  x<-reshape( data, idvar=idvar, varying=varying, direction="long", 
              times=times, v.names=v.names, timevar=within.tmp,
              new.row.names = paste0("case.",1:(dim(data)[1]*length(times))) )

  
  if( split==TRUE & length( grep(sep,times,fixed=TRUE))>0 ) { # split multiple treatments into two factors?
    split.treatments <- t(as.data.frame(strsplit(x[,within.tmp],sep,fixed=TRUE)))
    rownames(split.treatments)<-NULL
    split.treatments <- as.data.frame(split.treatments)
    if( length(within)==1) { 
      names(split.treatments) <- paste(within,1:length(split.treatments),sep="") 
    } else {
      if( length(within) == length(split.treatments)) {
        names(split.treatments) <- within 
      } else { stop( "length of 'within' is incorrect" )}
    }
    x <- x[,setdiff(names(x),within.tmp)] # delete collapsed treatment
    x <- cbind(x,split.treatments) # append split treatment
  } else { 
    x[,within.tmp]<- factor(x[,within.tmp])
    names(x)[grep(within.tmp,names(x))] <- within 
  }
  rownames(x) <- NULL
  names(x) <- make.unique(names(x))
  return(x)
  
}

