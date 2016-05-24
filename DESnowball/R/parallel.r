start.para <-
function(ncore, varlist, type=NULL) {
  if(ncore>1) {
      if(is.null(type)) cl <- makeCluster(ncore)
      else cl <- makeCluster(ncore,type=type)
  }
  else cl <- NULL
 clusterExport(cl, 
		c(varlist), 
		envir=parent.frame())
  .dump <- clusterEvalQ(cl, eval(require(DESnowball)))
  cl
}

stop.para <-
function(cl) {
  stopCluster(cl)
}
#
#snowball.initexpr <- function(){
#  expression({
#      require(DESnowball)
#  })
#}
