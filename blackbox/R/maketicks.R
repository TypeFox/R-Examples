maketicks <- function(var, labels=NULL, axis, maxticks) {
  ## we want pretty labels so we must determine xat from pretty xlabs values, not xlabs from pretty xat values...
  if (is.null(labels)) labels <- gridfn(var)
  if(islogscale(var)) {
    labels <- exp(labels)
    locstring <- " on a log scale"
  } else locstring <- ""
  ## at this point labels are in canonical scale
  if(var=="latt2Ns2") {
    Nbfactor <- blackbox.getOption("Nbfactor")
    labels <- labels*Nbfactor ## should be replaced by more general method for KrigSpace
  }
  labels <- niceLabels(labels, log=islogscale(var), axis=axis, maxticks=maxticks)
  if(islogscale(var)) {
    phantomat <- niceLabels(c(min(labels)/10, max(labels)*10), log=TRUE, lpos=c(1:9))
  } else phantomat <- NULL
  if(var=="latt2Ns2") {labels <- labels[labels>0];phantomat <- phantomat[phantomat>0]}
  ## as noted above we determine the at values from the nice labels values
  at <- labels
  if(var=="latt2Ns2") {
    at <- at/Nbfactor
    phantomat <- phantomat/Nbfactor
  } ## should be replaced by more general method for KrigSpace
  if(islogscale(var)) {at <- log(at);phantomat <- log(phantomat)}
  legend <- userunit(var, locstring, format="expression")
  return(list(labels=labels, at=at, legend=legend, phantomat=phantomat))
}
