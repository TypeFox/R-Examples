dpearson <-
function(x,params,moments,log=FALSE,...) {
  if (!missing(moments)) {
    if (!missing(params)) warning("'params' ignored, because 'moments' are given")
    params <- pearsonFitM(moments=moments)
  }
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  fname  <- paste("dpearson",suffix[params[[1]]+1],sep="")
  do.call(fname,args=list(x=x,params=params[-1],log=log,...))
}

ppearson <-
function(q,params,moments,lower.tail=TRUE,log.p=FALSE,...) {
  if (!missing(moments)) {
    if (!missing(params)) warning("'params' ignored, because 'moments' are given")
    params <- pearsonFitM(moments=moments)
  }
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  fname  <- paste("ppearson",suffix[params[[1]]+1],sep="")
  do.call(fname,args=list(q=q,params=params[-1],lower.tail=lower.tail,log.p=log.p,...))
}

qpearson <-
function(p,params,moments,lower.tail=TRUE,log.p=FALSE,...) {
  if (!missing(moments)) {
    if (!missing(params)) warning("'params' ignored, because 'moments' are given")
    params <- pearsonFitM(moments=moments)
  }
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  fname  <- paste("qpearson",suffix[params[[1]]+1],sep="")
  do.call(fname,args=list(p=p,params=params[-1],lower.tail=lower.tail,log.p=log.p,...))
}

rpearson <-
function(n,params,moments,...) {
  if (!missing(moments)) {
    if (!missing(params)) warning("'params' ignored, because 'moments' are given")
    params <- pearsonFitM(moments=moments)
  }
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  fname  <- paste("rpearson",suffix[params[[1]]+1],sep="")
  do.call(fname,args=list(n=n,params=params[-1],...))
}

pearsonMoments <-
function(params,moments) {
  if (!missing(moments)) {
    if (!missing(params)) warning("'params' ignored, because 'moments' are given")
    params <- pearsonFitM(moments=moments)
  }
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  fname  <- paste("pearson",suffix[params[[1]]+1],"moments",sep="")
  do.call(fname,args=params[-1])
}

