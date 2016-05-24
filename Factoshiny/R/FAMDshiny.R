FAMDshiny <-
function(X){
#   gassign("x", X)
   G <- .GlobalEnv
   assign("x", X, envir=G)
  nom=sys.calls()[[1]]
  nameJDD=nom[2]
#  gassign("nomData",nameJDD)
  assign("nomData",nameJDD, envir=G)
  if (!(inherits(X, "FAMDshiny") | inherits(X, "data.frame") | inherits(X, "FAMD"))){
    stop('df is not a dataframe, the results of the FAMDshiny function or a FAMD result')
  }
  if(is.data.frame(X)==TRUE){
    quanti=names(which(sapply(X,is.numeric)))
    quali=names(which(!(sapply(X,is.numeric))))
    if(length(quanti)==0 || length(quali)==0)
      stop('you data is not mixed')
  }
  a=shiny::runApp(system.file("FactoFAMDapp2", package="Factoshiny"),launch.browser = TRUE)
  return(invisible(a))
}



