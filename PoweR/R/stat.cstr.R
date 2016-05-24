stat.cstr <- function(stat.index,stat.pars=NULL) {

  if(getRversion() < "3.1.0") dontCheck <- identity

# We get the default (maximal) number of parameters
  Cstat.name <- paste("stat",stat.index,sep="")
  out <- .C(dontCheck(Cstat.name),0.0,0L,0.0,0L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,0L,0L,0.0,nbparams=1L)
  nbparams <- out$nbparams
  if (length(stat.pars) > nbparams) stop(paste("Length of 'stat.pars' should be at most",nbparams))
# We get the default values of the parameters (using the trick of putting the first value of *name to "1"
  out <- .C(dontCheck(Cstat.name),0.0,0L,0.0,0L,name=c("1",rep(" ",49)),1L,0.0,0L,0.0,0.0,0.0,0L,alter=0L,0L,params=rep(0.0,nbparams),nbparams=as.integer(nbparams))
  if (length(stat.pars) >= 1) stat.pars <- c(stat.pars,out$params[-(1:length(stat.pars))])
  if (is.null(stat.pars)) stat.pars <- out$params[1:nbparams]
  
  name <- out$name
  name <- gsub('\\','',gsub('$','',sub(' +$', '', paste(out$name,collapse="")),fixed=TRUE),fixed=TRUE)
  
  return(list(name=name,nbparams=nbparams,stat.pars=stat.pars,alter=if (out$alter == 0) c("0,1,2") else out$alter)) 
}
