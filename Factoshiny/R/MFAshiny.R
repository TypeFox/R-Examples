MFAshiny <-
function(X){
#    gassign("x", X)
    G <- .GlobalEnv
    assign("x", X, envir=G)

    if (!(inherits(X, "MFA") | inherits(X, "data.frame") | inherits(X, "MFAshiny"))){
      stop('df is not a dataframe or the result of the MFA function')
    }
    if(is.data.frame(X)==TRUE){
      nom=sys.calls()[[1]]
      nameJDD=nom[2]
#      gassign("nomData",nameJDD)
      assign("nomData",nameJDD, envir=G)
      if(dim(X)[2]<=2)
        stop('not enough variables in your dataset')
    a=shiny::runApp(system.file("FactoMFAapp", package="Factoshiny"),launch.browser = TRUE)
    return(invisible(a))
    }
    else{
      nom=as.character(X$call$call)[2] 
#      gassign("nomData",nom)
      assign("nomData",nom, envir=G)
a=shiny::runApp(system.file("FactoMFAapp2", package="Factoshiny"),launch.browser = TRUE)
      return(invisible(a))
    }
  }
