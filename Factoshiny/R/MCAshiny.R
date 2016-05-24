MCAshiny <-
  function(X){
#    gassign("x", X)
    G <- .GlobalEnv
    assign("x", X, envir=G)
    nom=sys.calls()[[1]]
    nameJDD=nom[2]
#    gassign("nomData", nameJDD)
    assign("nomData", nameJDD, envir=G)
    
    if((is.data.frame(X)==FALSE)&&(class(X)!="MCAshiny")&&!(inherits(X, "MCA")))
      stop('X is not a dataframe or a result from a MCA analysis')
    
   
    
    if((class(X)=="MCAshiny")||(is.data.frame(X)==TRUE)||(inherits(X, "MCA"))){
      ###
      if(is.data.frame(X)==TRUE){
        quali=names(which(!(sapply(X,is.numeric))))
      if(length(quali)<=2)
        stop('not enought qualitative variables in your dataset') 
      }
      ###
      a=shiny::runApp(system.file("FactoMCAapp2",package="Factoshiny"),launch.browser = TRUE)
      return(invisible(a))
    }
  }
