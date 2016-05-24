### Class definitions ###
## Types for representations
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("tsOrNULL", c("ts", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("listOrNULLOrnumeric", c("list", "numeric","NULL"))
setClassUnion("listOrNULLOrcharacter", c("list", "character","NULL"))
setClassUnion("numericOrNULLOrcharacter", c("numeric", "character","NULL"))
### Parameter Object class: x12Parameter ###
setClass(
    Class="x12Parameter", 
    representation=representation(			
        #period="numeric",
        series.span="numericOrNULLOrcharacter",
        series.modelspan="numericOrNULLOrcharacter",
        #series.type="characterOrNULL",
        #decimals="numeric",
        transform.function="character",
        transform.power="numericOrNULL",
        transform.adjust="characterOrNULL",
        regression.variables="characterOrNULL",
        regression.user="characterOrNULL",
        regression.file="characterOrNULL",
        regression.usertype="characterOrNULL",
        regression.centeruser="characterOrNULL",
        regression.start="numericOrNULLOrcharacter",
        regression.aictest="characterOrNULL",
        #outlier="logical",
        outlier.types="characterOrNULL",
        outlier.critical="listOrNULLOrnumeric",
        outlier.span="numericOrNULLOrcharacter",
        outlier.method="characterOrNULL",	
        identify="logical",
        identify.diff="numericOrNULL",
        identify.sdiff="numericOrNULL",
        identify.maxlag="numericOrNULL",	
        arima.model="numericOrNULL",
        arima.smodel="numericOrNULL",
        arima.ar="numericOrNULLOrcharacter",
        arima.ma="numericOrNULLOrcharacter",
        automdl="logical",
        automdl.acceptdefault="logical",
        automdl.balanced="logical",
        automdl.maxorder="numeric",
        automdl.maxdiff="numeric",
        forecast_years="numericOrNULL",
        backcast_years="numericOrNULL",
        forecast_conf="numeric",
        estimate="logical",
        estimate.outofsample="logical",
        check="logical",
        check.maxlag="numericOrNULL",
        slidingspans="logical",
        slidingspans.fixmdl="characterOrNULL",
        slidingspans.fixreg="characterOrNULL",
        slidingspans.length="numericOrNULL",
        slidingspans.numspans="numericOrNULL",
        slidingspans.outlier="characterOrNULL",
        slidingspans.additivesa="characterOrNULL",
        slidingspans.start="numericOrNULLOrcharacter",
        history="logical",
        history.estimates="characterOrNULL",
        history.fixmdl="logical",
        history.fixreg="characterOrNULL",
        history.outlier="characterOrNULL",
        history.sadjlags="numericOrNULL",
        history.trendlags="numericOrNULL",
        history.start="numericOrNULLOrcharacter",
        history.target="characterOrNULL",
        x11.sigmalim="numericOrNULL",
        x11.type="characterOrNULL",#vorher onlytd="logical"
        x11.sfshort="logical",
        x11.samode="characterOrNULL",
        x11.seasonalma="characterOrNULL",
        x11.trendma="numericOrNULL",
        x11.appendfcst="logical",
        x11.appendbcst="logical",
        x11.calendarsigma="characterOrNULL",
        x11.excludefcst="logical", 
        x11.final="character",
        x11regression="logical"
#	seats="logical",
#	seatsparameter="characterOrNULL"
    ),
    prototype=prototype(
        series.span=NULL,
        series.modelspan=NULL,
        #series.type=NULL,
        transform.function="auto",
        transform.power=NULL,
        transform.adjust=NULL,
        regression.variables=NULL,
        regression.user=NULL,
        regression.file=NULL,
        regression.usertype=NULL,
        regression.centeruser=NULL,
        regression.start=NULL,
        regression.aictest=NULL,
        #outlier=FALSE,
        outlier.types=NULL,
        outlier.critical=NULL,
        outlier.span=NULL,
        outlier.method=NULL,
        identify=FALSE,
        identify.diff=NULL,
        identify.sdiff=NULL,
        identify.maxlag=NULL,
        arima.model=NULL,
        arima.smodel=NULL,
        arima.ar=NULL,
        arima.ma=NULL,
        automdl=TRUE,
        automdl.acceptdefault=FALSE,
        automdl.balanced=TRUE,
        automdl.maxorder=c(3,2),
        automdl.maxdiff=c(1,1),
        forecast_years=1,
        backcast_years=NULL,
        forecast_conf=.95,
        estimate=FALSE,
        estimate.outofsample=TRUE,
        check=TRUE,
        check.maxlag=NULL,
        slidingspans=FALSE,
        slidingspans.fixmdl=NULL,
        slidingspans.fixreg=NULL,
        slidingspans.length=NULL,
        slidingspans.numspans=NULL,
        slidingspans.outlier=NULL,
        slidingspans.additivesa=NULL,
        slidingspans.start=NULL,
        history=FALSE,
        history.estimates=NULL,
        history.fixmdl=FALSE,
        history.fixreg=NULL,
        history.outlier=NULL,
        history.sadjlags=NULL,
        history.trendlags=NULL,
        history.start=NULL,
        history.target=NULL,
        x11.sigmalim=c(1.5,2.5),
        x11.type=NULL,
        x11.sfshort=FALSE,
        x11.samode=NULL,
        x11.seasonalma=NULL,
        x11.trendma=NULL,
        x11.appendfcst=TRUE,
        x11.appendbcst=FALSE,
        x11.calendarsigma=NULL,
        x11.excludefcst=FALSE,
        x11.final="user",
        x11regression=FALSE
#			seats=FALSE, 
#			seatsparameter=NULL
    ),
    validity=function(object) {
      return(TRUE)
    }
)
setClass(
    Class="spectrum", 
    representation=representation(
        frequency="numeric",
        spectrum="numeric"
    ),prototype=
        prototype(
            frequency=new("numeric"),
            spectrum=new("numeric")
        ),
    validity=function(object) {
      length(object@spectrum)==length(object@frequency)
    }
)
setClass(
    Class="fbcast",
    representation=representation(
        estimate="ts",
        lowerci="ts",
        upperci="ts"
    ),prototype=
        prototype(
            estimate=new("ts"),
            lowerci=new("ts"),
            upperci=new("ts")
        ),
    validity=function(object) {
      length(object@estimate)==length(object@lowerci)&&length(object@estimate)==length(object@upperci)
    }
)
setClass(
    Class="x12BaseInfo",
    representation=representation(
        x12path = "characterOrNULL",
        x13path = "characterOrNULL",
        use = "character",
        showWarnings = "logical" 
    ),prototype=
        prototype(
            x12path = NULL,
            x13path = NULL,
            use = "x12",
            showWarnings = FALSE
        ),
    validity=function(object) {
      (!is.null(object@x12path)||!is.null(object@x13path))&&object@use%in%c("x12","x13")
    }
)
setClass(Class="diagnostics",contains="list")
### Output Object class: x12Output ###
setClass(
    Class="x12Output", 
    representation=representation(			
        a1="ts",
        d10="ts",
        d11="ts",
        d12="ts",
        d13="ts",
        d16="ts",
        c17="ts",
        d9="ts",
        e2="ts",
        d8="ts",
        b1="ts",
        td="tsOrNULL",
        otl="tsOrNULL",
        sp0="spectrum",
        sp1="spectrum",
        sp2="spectrum",
        spr="spectrum",
        forecast="fbcast",
        backcast="fbcast",
        dg="list",
#      seats="logical",
        file="character",
        tblnames="character",
        Rtblnames="character"
    ),
    prototype=prototype(
        a1=new("ts"),
        d10=new("ts"),
        d11=new("ts"),
        d12=new("ts"),
        d13=new("ts"),
        d16=new("ts"),
        c17=new("ts"),
        d9=new("ts"),
        e2=new("ts"),
        d8=new("ts"),
        b1=new("ts"),
#		td=new("ts"),
#		otl=new("ts"),		
        sp0=new("spectrum"),
        sp1=new("spectrum"),
        sp2=new("spectrum"),
        spr=new("spectrum"),
        forecast=new("fbcast"),
        backcast=new("fbcast"),
        dg=new("diagnostics"),
#        seats=new("logical"),
        file=new("character"),
        tblnames=new("character"),
        Rtblnames=new("character")        
    ),
    validity=function(object) {
      return(TRUE)
    }
)

setClass(
    Class="x12Single", 
    representation=representation(			
        ts="ts",
        x12Parameter="x12Parameter",
        x12Output="x12Output",
        x12OldParameter="list",
        x12OldOutput="list",
        
        tsName="characterOrNULL",
        firstRun="logical"
    ),
    prototype=prototype(
        ts=new("ts"),
        x12Parameter=new("x12Parameter"),
        x12Output=new("x12Output"),
        x12OldParameter=new("list"),
        x12OldOutput=new("list"),
        tsName=NULL,
        firstRun=FALSE
    ),
    validity=function(object) {
      return(TRUE)
    }
)
setClass(Class="x12List",contains="list",validity=function(object){
      all(lapply(object,class),"x12Single")      
    })
setClass(
    Class="x12Batch", 
    representation=representation(			
        x12List="x12List",
        x12BaseInfo="x12BaseInfo"
    ),
    prototype=prototype(
        x12List=new("x12List"),
        x12BaseInfo=new("x12BaseInfo",use="x12",x12path="x12adummy")
    ),
    validity=function(object) {
      return(TRUE)
    }
)

setMethod(
    f='initialize',
    signature=signature(.Object = "x12Batch"),
    definition=function(.Object,tsList,tsName=NULL,x12BaseInfo=new("x12BaseInfo")) {
      for(i in 1:length(tsList)){
        if(class(tsList[[i]])=="x12Single")
          .Object@x12List[[i]] <- tsList[[i]]
        else{
          if(!is.null(tsName))
            .Object@x12List[[i]] <-new("x12Single",ts=tsList[[i]],tsName=tsName[i])
          else{
            .Object@x12List[[i]] <-new("x12Single",ts=tsList[[i]],tsName=paste("Series_",i,sep=""))
          }
        }
      }
      .Object@x12BaseInfo <- x12BaseInfo
      return(.Object)
    }
)
###Handling of x12path
setMethod(
    f='initialize',
    signature=signature(.Object = "x12BaseInfo"),
    definition=function(.Object,x12path=NULL,x13path=NULL,use=NULL,showWarnings=FALSE) {
      if(is.null(x12path)&&is.null(x13path)&&!existd("x12path")&&!existd("x13path"))
        stop("Please use the functions x12path() or x13path() to define the paths to the binaries.")
      if(is.null(x12path)&&existd("x12path")){
        if(file.exists(getd("x12path")))
          .Object@x12path <- getd("x12path")
        else
          stop("file specified in global variable x12path does not exist!\n")
      }
      if(is.null(x13path)&&existd("x13path")){
        if(file.exists(getd("x13path")))
          .Object@x13path <- getd("x13path")
        else
          stop("file specified in global variable x13path does not exist!\n")
      }
      if(!is.null(x12path)){
        if(file.exists(x12path)||x12path=="x12adummy")
          .Object@x12path <- x12path
        else 
          stop("file specified in argument x12path does not exist!\n")
      }
      if(!is.null(x13path)){
        if(file.exists(x13path))
          .Object@x13path <- x13path
        else
          stop("file specified in argument x13path does not exist!\n")
      }
      if(is.null(.Object@x12path)&&is.null(.Object@x13path)){
        stop("Please use the functions x12path() or x13path() to define the paths to the binaries.")
      }else if(!is.null(.Object@x12path)&&!is.null(.Object@x13path)){
        warning("x12path and x13path defined, x13 will be used!")
      }else if(!is.null(.Object@x13path)){
        use <- "x12"
        .Object@x12path <- .Object@x13path
        .Object@x13path <- NULL
      }else if(!is.null(.Object@x12path)){
        use <- "x12"
      }
      .Object@use <- use
      .Object@showWarnings <- showWarnings
      return(.Object)
    }
)
#Basic methods for x12Batch and x12Single
setMethod(
    f='dim',
    signature=signature(x = "x12Batch"),
    definition=function(x) {
          return(length(x@x12List))
    }
)
setMethod(
    f='length',
    signature=signature(x = "x12Batch"),
    definition=function(x) {
      return(length(x@x12List))
    }
)

setMethod(
    f='print',
    signature=signature(x = "x12Batch"),
    definition=function(x) {
      cat("A batch of time series of length ",length(x@x12List),".\n")
      for(i in 1:length(x@x12List))
        print(x@x12List[[i]])
    }
)

setMethod(
    f='print',
    signature=signature(x = "x12Single"),
    definition=function(x) {
      cat("Name: ",x@tsName,"\n")
      cat("processed with x12: ",x@firstRun,"\n")
      print(x@ts)
    }
)



setClass(Class="crossValidation",
		representation=representation(
				backcast="dfOrNULL",
				forecast="dfOrNULL"),
		prototype=prototype(
				backcast=NULL,
				forecast=NULL),
#		validity=function(object) {
#			length(object@spectrum)==length(object@frequency)
#		}
)
