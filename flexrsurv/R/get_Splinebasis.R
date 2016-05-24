# functions get_Splinebasis to get x spline basis
# functions get_TimeSplinebasis to get time spline basis

get_Splinebasis <- function(objterm,
                            data=parent.frame(),
                            specials="NPHNLL",
                            all.vars.func=all_specials_vars, 
                            unique=TRUE,
                            order=c("formula", "specials")){
  # get spline parameters of each NPHNLL terms
# input
#      terms : a term object 
# output  : list of "SplineBasis" objects

  order <- match.arg(order)
  
  indxvar <- attr(objterm, "specials")[specials]
  nvars <- length(unlist(indxvar))
  
  if(nvars==0){
    # no "specials" vars 
    return(NULL)
  }
  else{
    if(order=="specials"){
      oindxvar <- 1:nvars
    }
    else {
      oindxvar <- order(unlist(indxvar))
    }  
    var_list <- NULL
    Spline_list <- NULL

    for(is in specials){
      fun <- mget(is,
              mode = "function",
              envir = parent.frame(), inherits=TRUE,
              ifnotfound=list(NULL))[[1]]
      for( i in indxvar[[is]]){
        thecall <- match.call(fun, attr(objterm,"variables")[[i+1]])
        
        thevar <- thecall[["x"]]
        
        Knots <- eval(as.expression(thecall[["Knots"]]))
        
        if( !is.null(thecall[["Boundary.knots"]]) ){
          therange <- eval(as.expression(thecall[["Boundary.knots"]]))
        }
        else {
            # compute the range of the variable 
          therange <- eval(call("range", thevar), envir=data)
      }        
      
        thecall[["Spline"]] <- ifelse(is.null(thecall[["Spline"]]),
                                      eval(formals(fun)$Spline)[1],
                                      thecall[["Spline"]])
        if( is.null(thecall[["Spline"]])){
                                        # default is b-spline
        
          thespline <- MSplineBasis(knots=c(therange[1],
                                      eval(as.expression(thecall[["Knots"]])), 
                                      therange[2]),
                                    degree=ifelse(is.null(thecall[["Degree"]]),
                                      formals(fun)[["Degree"]],
                                      thecall[["Degree"]]), 
                                    keep.duplicates=FALSE)
        }
        else if( thecall[["Spline"]]== "tp-spline" ){
          thespline <- TPSplineBasis(knots=eval(as.expression(thecall[["Knots"]])), 
                                     degree=ifelse(is.null(thecall[["Degree"]]),
                                       formals(fun)[["Degree"]],
                                       thecall[["Degree"]]), 
                                     min=therange[1],
                                     max=therange[2],
                                     type="standard")
        }
        else if( thecall[["Spline"]]== "tpi-spline" ){
          thespline <- TPSplineBasis(knots=eval(as.expression(thecall[["Knots"]])), 
                                     degree=ifelse(is.null(thecall[["Degree"]]),
                                       formals(fun)[["Degree"]],
                                       thecall[["Degree"]]), 
                                     min=therange[1],
                                     max=therange[2],
                                     type="increasing")
        }
        else if( thecall[["Spline"]]== "b-spline" ){
          if (is.null(thecall[["Degree"]])) {thecall[["Degree"]]<-3}
          
          thespline <- MSplineBasis(knots=c(therange[1],
                                      eval(as.expression(thecall[["Knots"]])), 
                                      therange[2]),
                                    degree=ifelse(is.null(thecall[["Degree"]]),
                                      formals(fun)[["Degree"]],
                                      thecall[["Degree"]]),
                                    keep.duplicates=FALSE)
        }
        else { 
          stop("wrong type of spline specification", attr(objterm,"variables")[[i+1]])
        }
        var_list <- c( var_list, thevar) 
        Spline_list <- c( Spline_list, thespline)
      }
    }
    names(Spline_list) <- var_list
    return(Spline_list[oindxvar])
  }
  
  
}



get_TimeSplinebasis <- function(objterm,
                                data=parent.frame(),
                                specials="NPHNLL",
                                all.vars.func=all_specials_vars, 
                                unique=TRUE,
                                order=c("formula", "specials")){
  # get spline parameters of each NPHNLL terms
# input
#      terms : a term object 
# output  : list of "SplineBasis" objects

  order <- match.arg(order)
  
  indxvar <- attr(objterm, "specials")[specials]
  nvars <- length(unlist(indxvar))
  
  if(nvars==0){
    # no "specials" vars 
    return(NULL)
  }
  else{
    if(order=="specials"){
      oindxvar <- 1:nvars
    }
    else {
      oindxvar <- order(unlist(indxvar))
    }  
    var_list <- NULL
    Spline_list <- NULL

    for(is in specials){
      fun <- mget(is,
              mode = "function",
              envir = parent.frame(), inherits=TRUE,
              ifnotfound=list(NULL))[[1]]
      for( i in indxvar[[is]]){
        thecall <- match.call(fun, attr(objterm,"variables")[[i+1]])
     
        thevar <- thecall[["timevar"]]
        
        Knots <- eval(as.expression(thecall[["Knots.t"]]))
        
        if( !is.null(thecall[["Boundary.knots.t"]]) ){
          therange <- eval(as.expression(thecall[["Boundary.knots.t"]]))
        }
        else {
                                        # compute the range of the variable 
          therange <- eval(call("range", thevar), envir=data)
        }        
        
        thecall[["Spline"]] <- ifelse(is.null(thecall[["Spline"]]),
                                      eval(formals(fun)$Spline)[1],
                                      thecall[["Spline"]])
        if( is.null(thecall[["Spline"]])){
                                        # default is b-spline
          
          thespline <- MSplineBasis(knots=c(therange[1],
                                      eval(as.expression(thecall[["Knots.t"]])), 
                                      therange[2]),
                                    degree=ifelse(is.null(thecall[["Degree.t"]]),
                                      formals(fun)[["Degree.t"]],
                                      thecall[["Degree.t"]]), 
                                    keep.duplicates=FALSE)
        }
        else if( thecall[["Spline"]]== "tp-spline" ){
          thespline <- TPSplineBasis(knots=eval(as.expression(thecall[["Knots.t"]])), 
                                     degree=ifelse(is.null(thecall[["Degree.t"]]),
                                       formals(fun)[["Degree.t"]],
                                       thecall[["Degree.t"]]), 
                                     min=therange[1],
                                     max=therange[2],
                                     type="standard")
        }
        else if( thecall[["Spline"]]== "tpi-spline" ){
          thespline <- TPSplineBasis(knots=eval(as.expression(thecall[["Knots.t"]])), 
                                     degree=ifelse(is.null(thecall[["Degree.t"]]),
                                       formals(fun)[["Degree.t"]],
                                       thecall[["Degree.t"]]), 
                                     min=therange[1],
                                     max=therange[2],
                                     type="standard")
        }
        else if( thecall[["Spline"]]== "b-spline" ){
          if (is.null(thecall[["Degree.t"]])) {thecall[["Degree.t"]]<-3}
          
          thespline <- MSplineBasis(knots=c(therange[1],
                                      eval(as.expression(thecall[["Knots.t"]])), 
                                      therange[2]),
                                    degree=ifelse(is.null(thecall[["Degree.t"]]),
                                      formals(fun)[["Degree.t"]],
                                      thecall[["Degree.t"]]),
                                    keep.duplicates=FALSE)
        }
        else { 
          stop("wrong type of spline specification", attr(objterm,"variables")[[i+1]])
        }
        var_list <- c( var_list, thevar) 
        Spline_list <- c( Spline_list, thespline)
      }
    }
    names(Spline_list) <- var_list
    return(Spline_list[oindxvar])
  }
  
  
}
