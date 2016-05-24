# individuals 14_07_24

individuals <- function(eloobject, from=eloobject$misc["maxDate"], to=NULL, outp=c("N", "IDs", "CV")) {
  # outp:
  # N - (mean) number of individuals
  # IDs - IDs that were present on at least one of the day(s)
  # CV - coefficient of variation of N
  
  # some checks and prelims
  outp <- match.arg(outp)
  
  if(!is.null(to)) { 
    if(to==from) to <- NULL 
    #  if(as.Date(to)<as.Date(from)) stop("the 'to' date lies before the starting ('from') date")
  }
  
  
  # create vector with all dates (according to date range in eloobject)
  DR <- seq(from=as.Date(eloobject$misc["minDate"]), to=as.Date(eloobject$misc["maxDate"]), by="day")
  
  # presence matrix
  pmat <- eloobject$pmat
  
  # if no 'to' date is given, i.e. only one day (also the case if 'to' is the same as 'from')
  if(is.null(to)) { 
    l <- which(DR==as.Date(from))
    if(outp=="N")   { res <- sum(pmat[l, ]) }
    if(outp=="CV")  { res <- NA }
    if(outp=="IDs") { res <- names(pmat[l, which(pmat[l, ]==1)]) }
  } 
  
  # date range is given, i.e. 'to' is other than NULL (and other than 'from')
  if(!is.null(to)) { 
    l <- which(DR %in% seq(from=as.Date(from), to=as.Date(to), by="day"))
    if(outp=="N")   { res <- mean(rowSums(pmat[l, ])) }
    if(outp=="CV")  { res <- sd(rowSums(pmat[l, ]))/mean(rowSums(pmat[l, ])) }
    if(outp=="IDs") { res <- names(which(colSums(pmat[l, ]) >= 1)) }    
  }
  
  
  return(res)
  
}

