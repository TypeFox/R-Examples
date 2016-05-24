law.cstr <- function(law.index,law.pars=NULL) {

  if(getRversion() < "3.1.0") dontCheck <- identity

# We get the default (maximal) number of parameters and the default values of the parameters
# We use the trick to put the first element of *name to "1" to also retrieve the default values
  Claw.name <- paste("law",law.index,sep="")
  out <- .C(dontCheck(Claw.name),0L,0.0,name=c("1",rep(" ",49)),1L,params=rep(0.0,4),nbparams=0L,0L)
  nbparams <- out$nbparams
  if (length(law.pars) > nbparams) stop(paste("Length of 'law.pars' should be at most",nbparams))
  if (length(law.pars) >= 1) law.pars <- c(law.pars,out$params[-(1:length(law.pars))])[1:nbparams]
  if (is.null(law.pars)) law.pars <- out$params[1:nbparams]
  
  name <- out$name
  name <- gsub('\\','',gsub('$','',sub(' +$', '', paste(name,collapse="")),fixed=TRUE),fixed=TRUE)
  if (length(grep("[()]",name)) == 1) {
    split1 <- unlist(strsplit(name, "(", fixed = TRUE))
    law.name <- split1[1]
    law.args <- unlist(strsplit(unlist(strsplit(split1[2],")", fixed = TRUE)),",", fixed = TRUE))
    name <- paste(law.name,"(",paste(law.args,law.pars,sep="=",collapse=","),")",sep="")
  } else {
    name <- paste(name,"(",paste(law.pars,sep="",collapse=","),")",sep="")
  }
  
  return(list(name=name,nbparams=nbparams,law.pars=law.pars)) 
}
