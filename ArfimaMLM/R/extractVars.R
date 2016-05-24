#' @export
extractVars <-
function(formula,ecmformula=NULL){
  
  # select original varnames for transformation
  vars <- all.vars(formula)
  vars.ecm <- all.vars(ecmformula)
  ydif <- sub(".ydif", "", vars[grep(".ydif", vars, fixed=T)], fixed=T)
  xdif <- unique(sub(".xdif", "", vars[grep(".xdif", vars, fixed=T)], fixed=T))
  xfd <- sub(".fd", "", vars[grep(".fd", vars, fixed=T)], fixed=T)
  ecm <- sub(".mean", "", vars.ecm[grep(".mean", vars.ecm, fixed=T)], fixed=T)
  
  # generate varlists for analyses
  varlist.mean <- unique(c(ydif,xdif,xfd,ecm))
  varlist.fd <- unique(c(ydif,xfd))
  
  # output
  out <- list(mean=varlist.mean,fd=varlist.fd,xdif=xdif,ydif=ydif)
  return(out)
}