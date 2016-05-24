setGeneric("x12",
    function(object, x12Parameter=new("x12Parameter"),
        x12BaseInfo=new("x12BaseInfo"),...) { standardGeneric("x12")} )
setMethod(
    f='x12',
    signature=signature(object = "ts"),
    definition=function(object, x12Parameter,x12BaseInfo) {
      Par <- slotNames(x12Parameter)
      pp <- vector()
      for(p in Par){
        pp <- c(pp,(paste(p,"=x12Parameter@",p,sep="")))
      }
      Par <- slotNames(x12BaseInfo)
       for(p in Par){
         pp <- c(pp,(paste(p,"=x12BaseInfo@",p,sep="")))
       }
       if(!is.null(getOption("x12.delete"))){
         if(getOption("x12.delete"))
           keep_x12out <- paste("keep_x12out=FALSE")
         else
           keep_x12out <- paste("keep_x12out=TRUE")
       }else
         keep_x12out <- paste("keep_x12out=TRUE")
      pp <- paste("out <- x12work(tso=object,",paste(pp,collapse=","),",tblnames=\"otl\",Rtblnames=\"regressor\",",keep_x12out,")",sep="")
      eval(parse(text=pp))
      classout <- new("x12Output")
      Par <- slotNames(classout)
       for(p in Par){
         if(class(slot(classout,p))=="spectrum"){
           if(p%in%names(out)){
             slot(classout,p)@frequency <- out[[p]]$frequency
             slot(classout,p)@spectrum <- out[[p]]$spectrum
           }
         }else if(class(slot(classout,p))=="fbcast"){
           if(p%in%names(out)){
             slot(classout,p)@estimate <- out[[p]][["estimate"]]
             slot(classout,p)@lowerci <- out[[p]][["lowerci"]]
             slot(classout,p)@upperci <- out[[p]][["upperci"]]
           }
         }else
           slot(classout,p)<-out[[p]] 
       }
      return(classout)
    }
)
setMethod(
    f='x12',
    signature=signature(object = "x12Single"),
    definition=function(object,x12BaseInfo=new("x12BaseInfo"),forceRun=FALSE) {
      if(length(object@x12OldParameter)>0)
        TF <- !identical(object@x12Parameter,object@x12OldParameter[[length(object@x12OldParameter)]])
      else
        TF <- TRUE
      if(!object@firstRun||forceRun||TF){
        x12Parameter <- object@x12Parameter  
        if(object@firstRun){
          object@x12OldParameter[[length(object@x12OldParameter)+1]] <- object@x12Parameter
          object@x12OldOutput[[length(object@x12OldOutput)+1]] <- object@x12Output
        }
        object@firstRun <- TRUE
        Par <- slotNames(x12Parameter)
        pp <- vector()
        for(p in Par){
          pp <- c(pp,(paste(p,"=x12Parameter@",p,sep="")))
        }
        Par <- slotNames(x12BaseInfo)
        for(p in Par){
          pp <- c(pp,(paste(p,"=x12BaseInfo@",p,sep="")))
        }
        
        if(!is.null(object@tsName))
          pp <- c(pp, paste("file=\"",object@tsName,"\"",sep=""))
        if(!is.null(getOption("x12.delete"))){
          if(getOption("x12.delete"))
            keep_x12out <- paste("keep_x12out=FALSE")
          else
            keep_x12out <- paste("keep_x12out=TRUE")
        }else
          keep_x12out <- paste("keep_x12out=TRUE")
        pp <- paste("out <- x12work(tso=object@ts,",paste(pp,collapse=","),",tblnames=\"otl\",Rtblnames=\"regressor\",",keep_x12out,")",sep="")
        eval(parse(text=pp))
        classout <- new("x12Output")
        Par <- slotNames(classout)
        for(p in Par){
          if(class(slot(classout,p))=="spectrum"){
            if(p%in%names(out)){
              slot(classout,p)@frequency <- out[[p]]$frequency
              slot(classout,p)@spectrum <- out[[p]]$spectrum
            }
          }else if(class(slot(classout,p))=="fbcast"){
            if(p%in%names(out)){
              slot(classout,p)@estimate <- out[[p]][["estimate"]]
              slot(classout,p)@lowerci <- out[[p]][["lowerci"]]
              slot(classout,p)@upperci <- out[[p]][["upperci"]]
            }
          }else
            slot(classout,p)<-out[[p]] 
        }
        object@x12Output <- classout
        
      }
      return(object)
    }
)
setMethod(
    f='x12',
    signature=signature(object = "x12Batch"),
    definition=function(object,forceRun=FALSE) {
      starting.time <- Sys.time()
      if(existd("x12path"))
        object@x12BaseInfo@x12path <- getd("x12path")
      else
        stop("Please enter an x12path")
      ## Parallelization implemented after the pattern used in the survey package by Thomas Lumley.
      if (is.null(getOption("x12.parallel")) | !require("parallel", quietly=TRUE)){
        tmpList <- lapply(object@x12List,function(x)try(x12(x,x12BaseInfo=object@x12BaseInfo,forceRun=forceRun)))
      }else{
        tmpList <- mclapply(object@x12List,function(x)try(x12(x,x12BaseInfo=object@x12BaseInfo,forceRun=forceRun)),mc.cores=getOption("x12.parallel"))
      }
      for(i in 1:length(tmpList))
        object@x12List[[i]] <- tmpList[[i]] 
      print(Sys.time()-starting.time)
      return(object)
    }
)
