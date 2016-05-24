# Generic function setP, set Parameter
setGeneric("setP",
    function(object, ...) { standardGeneric("setP")} )
# Generic function getP, get Parameter
setGeneric("getP",
    function(object, ...) { standardGeneric("getP")} )
#Methods for signature x12Parameter
setMethod(
    f='setP',
    signature=signature(object = "x12Parameter"),
    definition=function(object, listP) {
      paras <- c(
			  #"period",
			  "series.span",
			  "series.modelspan",
			  #"series.type",
			  #"decimals",
			  "transform.function",
			  "transform.power",
			  "transform.adjust",
			  "regression.variables",
			  "regression.user",
			  "regression.file",
			  "regression.usertype",
			  "regression.centeruser",
			  "regression.start",
			  "regression.aictest",
			  #"outlier",
			  "outlier.types",
			  "outlier.critical",
			  "outlier.span",
			  "outlier.method",
			  "identify",
			  "identify.diff",
			  "identify.sdiff",
			  "identify.maxlag",
			  "arima.model",
			  "arima.smodel",
			  "arima.ar",
			  "arima.ma",
			  "automdl",
			  "automdl.acceptdefault",
			  "automdl.balanced",
			  "automdl.maxorder",
			  "automdl.maxdiff",
			  "forecast_years",
			  "backcast_years",
			  "forecast_conf",
			  "estimate",
			  "estimate.outofsample",
			  "check",
			  "check.maxlag",
			  "slidingspans",
			  "slidingspans.fixmdl",
			  "slidingspans.fixreg",
			  "slidingspans.length",
			  "slidingspans.numspans",
			  "slidingspans.outlier",
			  "slidingspans.additivesa",
			  "slidingspans.start",
			  "history",
			  "history.estimates",
			  "history.fixmdl",
			  "history.fixreg",
			  "history.outlier",
			  "history.sadjlags",
			  "history.trendlags",
			  "history.start",
			  "history.target",
			  "x11.sigmalim",
			  "x11.type",
			  "x11.sfshort",
			  "x11.samode",
			  "x11.seasonalma",
			  "x11.trendma",
			  "x11.appendfcst",
			  "x11.appendbcst",
			  "x11.calendarsigma",
			  "x11.excludefcst",
			  "x11.final",
			  "x11regression"
			  #"tblnames",
			  #"Rtblnames",
			  #"seats",
			  #"seatsparameter"
				)
      mn <- names(listP)%in%paras
      if(any(!mn)){
        warning("The following parameters could not be matched: ",paste(names(listP)[!mn],collapse=" , "))
      }
      mn <- names(listP)[mn]
      for(nam in mn){
        slot(object,nam) <- listP[[nam]]
      }
      return(object)
    }
)

setMethod(
    f='getP',
    signature=signature(object = "x12Parameter"),
    definition=function(object, whichP) {
      paras <- c(
			  #"period",
			  "series.span",
			  "series.modelspan",
			  #"series.type",
			  #"decimals",
			  "transform.function",
			  "transform.power",
			  "transform.adjust",
			  "regression.variables",
			  "regression.user",
			  "regression.file",
			  "regression.usertype",
			  "regression.centeruser",
			  "regression.start",
			  "regression.aictest",
			  #"outlier",
			  "outlier.types",
			  "outlier.critical",
			  "outlier.span",
			  "outlier.method",
			  "identify",
			  "identify.diff",
			  "identify.sdiff",
			  "identify.maxlag",
			  "arima.model",
			  "arima.smodel",
			  "arima.ar",
			  "arima.ma",
			  "automdl",
			  "automdl.acceptdefault",
			  "automdl.balanced",
			  "automdl.maxorder",
			  "automdl.maxdiff",
			  "forecast_years",
			  "backcast_years",
			  "forecast_conf",
			  "estimate",
			  "estimate.outofsample",
			  "check",
			  "check.maxlag",
			  "slidingspans",
			  "slidingspans.fixmdl",
			  "slidingspans.fixreg",
			  "slidingspans.length",
			  "slidingspans.numspans",
			  "slidingspans.outlier",
			  "slidingspans.additivesa",
			  "slidingspans.start",
			  "history",
			  "history.estimates",
			  "history.fixmdl",
			  "history.fixreg",
			  "history.outlier",
			  "history.sadjlags",
			  "history.trendlags",
			  "history.start",
			  "history.target",
			  "x11.sigmalim",
			  "x11.type",
			  "x11.sfshort",
			  "x11.samode",
			  "x11.seasonalma",
			  "x11.trendma",
			  "x11.appendfcst",
			  "x11.appendbcst",
			  "x11.calendarsigma",
			  "x11.excludefcst",
			  "x11.final",
			  "x11regression"
			  #"tblnames",
			  #"Rtblnames",
			  #"seats",
			  #"seatsparameter"
			  )
      
      mn <- whichP%in%paras
      if(any(!mn)){
        warning("The following parameters could not be matched: ",paste(whichP[!mn],collapse=" , "))
      }
      mn <- whichP[mn]
      ret <- list()
      for(nam in mn){
        ret[[nam]] <- slot(object,nam)
      }
      return(ret)
    }
)
#Methods for signature x12Single
setMethod(
    f='getP',
    signature=signature(object = "x12Single"),definition=function(object, whichP) {
      getP(object@x12Parameter,whichP=whichP)
    })
setMethod(
    f='setP',
    signature=signature(object = "x12Single"),definition=function(object, listP) {
      object@x12Parameter <- setP(object@x12Parameter,listP=listP)
      return(object)
    })
#Methods for signature x12Batch
setMethod(
    f='getP',
    signature=signature(object = "x12Batch"),definition=function(object, whichP,index=NULL) {
      ret <- list()
      if(is.null(index)){##changing all
        cat("The parameters for all objects are shown.\n")
        for(i in 1:length(object@x12List)){
          ret[[length(ret)+1]] <- getP(object@x12List[[i]],whichP=whichP)
        } 
      }else{
        if(is.integer(index)){
          if(min(index)>0&max(index)<=length(object@x12List)){
            for(i in index){
              ret[[length(ret)+1]] <- getP(object@x12List[[i]],whichP=whichP)
            }
          }else
            stop("argument index is out of bounds!\n")
        }else if(is.character(index)){
          namTS <- vector()
          for(i in 1:length(object@x12List)){
            namTS <- c(namTS,object@x12List[[i]]@tsName)
          }
          if(all(index%in%namTS)){
            for(nam in index){
              ind <- which(nam==namTS)
              ret[[length(ret)+1]] <- getP(object@x12List[[ind]],whichP=whichP)
            }
          }else
            stop("argument index contained names not found in the series names!\n")
          
        }else
          stop("argument index must be either integer or character!\n")
        
      }
      return(ret)
    })

setMethod(
    f='setP',
    signature=signature(object = "x12Batch"),definition=function(object, listP,index=NULL) {
      if(is.null(index)){##changing all
        cat("The parameters for all objects are changed.\n")
        for(i in 1:length(object@x12List)){
          object@x12List[[i]] <- setP(object@x12List[[i]],listP=listP)
        } 
      }else{
        if(is.numeric(index)){
          if(min(index)>0&max(index)<=length(object@x12List)){
            for(i in index){
              object@x12List[[i]] <- setP(object@x12List[[i]],listP=listP)
            }
          }else
            stop("argument index is out of bounds!\n")
        }else if(is.character(index)){
          namTS <- vector()
          for(i in 1:length(object@x12List)){
            namTS <- c(namTS,object@x12List[[i]]@tsName)
          }
          if(all(index%in%namTS)){
            for(nam in index){
              ind <- which(nam==namTS)
              object@x12List[[ind]] <- setP(object@x12List[[ind]],listP=listP)
            }
          }else
            stop("argument index contained names not found in the series names!\n")
          
        }else
          stop("argument index must be either integer or character!\n")

      }
      
      return(object)
    })
#Goto previous parameter setting and output
# Generic function prev, cleanArchive
setGeneric("prev",
    function(object, ...) { standardGeneric("prev")} )
setMethod(
    f='prev',
    signature=signature(object = "x12Single"),definition=function(object,n=NULL) {
      if(is.null(n))
        ind <- length(object@x12OldParameter)
      else if(n%in%c(1:length(object@x12OldParameter)))
        ind <- n
      else
        stop("Please provide an index corresponding to a previous run. (see summary with oldOutput>0)")
        
      object@x12Output <- object@x12OldOutput[[ind]]
      object@x12Parameter <- object@x12OldParameter[[ind]]
      oldout <- list()
      oldpar <- list()
      for(i in 1:length(object@x12OldParameter)){
        if(i!=ind){
          oldout[[length(oldout)+1]] <- object@x12OldOutput[[i]]
          oldpar[[length(oldpar)+1]] <- object@x12OldParameter[[i]]
        }
      }
      object@x12OldOutput <- oldout
      object@x12OldParameter <- oldpar
      return(object)
    })
setMethod(
    f='prev',
    signature=signature(object = "x12Batch"),definition=function(object,index=NULL,n=NULL) {
      if(is.null(index)){##changing all
        cat("All current parameters and outputs are replaced by the previous ones.\n")
        for(i in 1:length(object@x12List)){
          object@x12List[[i]] <- prev(object@x12List[[i]],n=n) 
        } 
      }else{
        if(is.numeric(index)){
          if(min(index)>0&max(index)<=length(object@x12List)){
            for(i in index){
              object@x12List[[i]] <- prev(object@x12List[[i]],n=n)
            }
          }else
            stop("argument index is out of bounds!\n")
        }else if(is.character(index)){
          namTS <- vector()
          for(i in 1:length(object@x12List)){
            namTS <- c(namTS,object@x12List[[i]]@tsName)
          }
          if(all(index%in%namTS)){
            for(nam in index){
              ind <- which(nam==namTS)
              object@x12List[[ind]] <- prev(object@x12List[[ind]],n=n)
            }
          }else
            stop("argument index contained names not found in the series names!\n")
          
        }else
          stop("argument index must be either integer or character!\n")
      }
      return(object)
    })
setGeneric("cleanArchive",
    function(object, ...) { standardGeneric("cleanArchive")} )
setGeneric("cleanHistory",
    function(object, ...) {
      .Deprecated("cleanArchive")
      cleanArchive(object,...)
    } )
setMethod(
    f='cleanArchive',
    signature=signature(object = "x12Single"),definition=function(object) {
       object@x12OldParameter <- object@x12OldOutput <- list()
      return(object)
    })
setMethod(
    f='cleanArchive',
    signature=signature(object = "x12Batch"),definition=function(object,index=NULL) {
      if(is.null(index)){##changing all
        cat("All previous parameters and outputs are deleted.\n")
        for(i in 1:length(object@x12List)){
          object@x12List[[i]] <- cleanArchive(object@x12List[[i]]) 
        } 
      }else{
        if(is.numeric(index)){
          if(min(index)>0&max(index)<=length(object@x12List)){
            for(i in index){
              object@x12List[[i]] <- cleanArchive(object@x12List[[i]])
            }
          }else
            stop("argument index is out of bounds!\n")
        }else if(is.character(index)){
          namTS <- vector()
          for(i in 1:length(object@x12List)){
            namTS <- c(namTS,object@x12List[[i]]@tsName)
          }
          if(all(index%in%namTS)){
            for(nam in index){
              ind <- which(nam==namTS)
              object@x12List[[ind]] <- cleanArchive(object@x12List[[ind]])
            }
          }else
            stop("argument index contained names not found in the series names!\n")
          
        }else
          stop("argument index must be either integer or character!\n")
      }
      return(object)
    })

####SAVE
setGeneric("saveP",
    function(object, file="x12Parameter.RData") { standardGeneric("saveP")} )
setGeneric("loadP",
    function(object, file) { standardGeneric("loadP")} )

setMethod(
    f='saveP',
    signature=signature(object = "x12Parameter"),
    definition=function(object,file) {
      save(object,file=file)
    }
)
setMethod(
    f='saveP',
    signature=signature(object = "x12Single"),
    definition=function(object,file) {
      out=object@x12Parameter
      save(out,file=file)
    }
)
setMethod(
    f='saveP',
    signature=signature(object = "x12Batch"),
    definition=function(object,file) {
      x12ParList <- list()
      for(i in 1:length(object@x12List)){
        x12ParList[[object@x12List[[i]]@tsName]] <- object@x12List[[i]]@x12Parameter
      }
      save(x12ParList,file=file)
    }
)
setMethod(
    f='loadP',
    signature=signature(object = "x12Parameter"),
    definition=function(object,file) {
      par <- get(load(file=file))
      if("x12Parameter"!=class(par))
        stop("no parameter settings found in the file!\n")
      return(par)
    }
)
setMethod(
    f='loadP',
    signature=signature(object = "x12Single"),
    definition=function(object,file) {
      par <- get(load(file=file))
      if("x12Parameter"!=class(par))
        stop("no parameter settings found in the file!\n")
      object@x12Parameter <- par
      return(object)
    }
)
setMethod(
    f='loadP',
    signature=signature(object = "x12Batch"),
    definition=function(object,file) {
      parList <- get(load(file=file))
      if(class(parList)=="x12Parameter"){
        warning("All Parameters will be overwritten with one loaded parameter configuration")
        for(i in 1:length(object@x12List)){
          object@x12List[[i]]@x12Parameter <- parList 
        }
      }else{
        if(length(parList)!=length(object@x12List))
          stop("loaded Parameter list does not fit to the x12Batch object \n")
        for(i in 1:length(parList)){
          if(class(parList[[i]])!="x12Parameter")
            stop("The file does not contain a list of x12Parameter objects!")
          object@x12List[[i]]@x12Parameter <- parList[[i]]
        }
      }
      return(object)
    }
)