##################################################################################

##' Class "NpdeData" representing the structure of the longitudinal data
##' 
##' A longitudinal data structure
##' 
##' @name NpdeData-class
##' @aliases NpdeData NpdeData-class show,NpdeData-method print,NpdeData-method 
##' summary,NpdeData-method
##' [,NpdeData-method [<-,NpdeData-method
##' @docType class
##' @section Objects from the Class: NpdeData objects are typically created by calls to \code{\link{npdeData}} and contain the following slots:
##' 
##' \describe{
##' \item{name.data}{character string giving the name of the dataset}
##' \item{name.group}{character string giving the name of the grouping term (ID)}
##' \item{name.predictor}{character string giving the name of the predictor (X)}
##' \item{name.response}{character string giving the name of the response (Y)}
##' \item{name.cens}{character string giving the name of the censoring indicator}
##' \item{name.mdv}{character string giving the name of the missing data indicator}
##' \item{name.covariates}{vector of character string giving the name(s) of the covariates}
##' \item{name.ipred}{character string giving the name of the individual predictions}
##' \item{units}{(optional) a list with the units for X, Y, and covariates}
##' \item{data}{a dataframe containing the data}
##' \item{N}{number of subjects}
##' \item{ntot.obs}{total number of non-missing observations}
##' \item{nind.obs}{vector of size N giving the number of non-missing observations for each subject}
##' \item{ind}{index of non-missing observations}
##' \item{icens}{index of censored observations (non-missing)}
##' \item{not.miss}{a vector of boolean indicating for each observation whether it is missing (FALSE) or available (TRUE)}
##' \item{loq}{the censoring value}
##' }
##' @section Methods:
##' \describe{
##'   \item{npdeData(name.data):}{Create a new \code{\linkS4class{NpdeData}} object from dataset name.data}
##'   \item{print(npde.data):}{Prints a summary of object npde.data}
##'   \item{show(npde.data):}{Prints a short summary of object npde.data}
##'   \item{showall(npde.data):}{Prints a detailed summary of object npde.data}
##'   \item{plot(npde.data):}{Plots the data in npde.data. More details can be found in \code{\link{plot.NpdeData}}}
##'   \item{summary(npde.data):}{Returns a summary of object npde.data in list format}
##'   \item{set.plotoptions(npde.data):}{Sets options for graphs of npde.data (internal method used in plots)}
##' }
##' @seealso \code{\link{npde}}, \code{\link{autonpde}}, \code{\link{plot.NpdeData}}
##' @keywords classes
##' @examples
##' 
##' methods(class="NpdeData")
##' 
##' showClass("NpdeData")
##' 
##' @exportClass NpdeData

setClass(
  Class="NpdeData",
  representation=representation(
    name.data="character",	# name of dataset
    name.group="character",	# name of column with ID
    name.predictor="character",# name of column(s) with predictors 
    name.response="character",	# name of column with response
    name.cens="character",	# name of column with censoring information
    name.miss="character",	# name of column indicating missing data (not censored)
    name.covariates="character",# name of column(s) with covariates (can be used for graphs)
    name.ipred="character",	# name of column indicating individual predictions (if available in the dataset)
    units="list",		# units (list with components for x, y, and cov), used in plots
    data="data.frame",          # the data: data frame with columns name.group (subject id), index (id renamed to 1:N), name.predictors (predictors), name.response (possibly transformed during fit), cens (1 if data is censored; the corresponding response is the censoring value), mdv (1 if the data is considered as missing), name.covariates (covariates, binary covariates are modified to 0/1), name.ipred (individual predictions, if available)
    ind="numeric",		# index of the non-missing observations # ECO TODO remove ?
    icens="numeric",		# index of the censored observations (non-missing)
    not.miss="logical",		# vector of logical, TRUE if present (=not missing), FALSE for missing data
    N="numeric",		# number of subjects
    ntot.obs="numeric",		# total number of observations
    nind.obs="numeric",		# number of observations for each subject
    loq="numeric"		# LOQ
    ),
  validity=function(object){
#    cat ("--- Checking NpdeData object ---\n")
    if (length(object@name.data)==0) {
      stop ("[ NpdeData : validation ] Please provide a name for the data (dataset or datafile on disk).")
    }
#     if(object@N>0) {
#       N<-object@N
#       if(N!=length(unique(object@id)) | N!=length(unique(object@index)) | N!=length(object@nind.obs)) {
#         cat("Size mismatch: id, index and/or nind.obs no longer correspond to the number of subjects.\n")
#       }
#       nobs<-object@ntot.obs
#       if(nobs!=dim(object@data)[1]) {
#         cat("Check length of predictor and response.\n")
#       }
#     }
    return(TRUE)
  }
)

##' Class "NpdeSimData" representing the structure of the longitudinal data
##' 
##' A longitudinal data structure, with simulated data
##' 
##' @name NpdeSimData-class
##' @aliases NpdeSimData NpdeSimData-class show,NpdeSimData-method [,NpdeSimData-method [<-,NpdeSimData-method
##' @docType class
##' @section Objects from the Class: NpdeSimData objects are created by associating an NpdeData object with matching simulated data, and they contain the following slots.
##' 
##' \describe{
##' \item{name.simdata}{character string giving the name of the dataset}
##' \item{nrep}{number of replications)}
##' \item{datsim}{a dataframe containing the simulated data,  with columns: idsim (subject id), irsim (replication index), xsim (simulated x), ysim (simulated response). After a call to \code{\link{npde}} or \code{\link{autonpde}}, an additional column ydsim (decorrelated replicated data) will be added.}
##' }
##' @section Methods:
##' \describe{
##'   \item{print(npde.simdata):}{Prints a summary of object npde.simdata}
##'   \item{show(npde.simdata):}{Prints a short summary of object npde.simdata}
##'   \item{showall(npde.simdata):}{Prints a detailed summary of object npde.simdata}
##' }
##' @seealso \code{\link{npde}}, \code{\link{autonpde}}
##' @keywords classes
##' @examples
##' 
##' showClass("NpdeSimData")
##' 
##' @export

setClass(
  Class="NpdeSimData",
  representation=representation(
    name.simdata="character",	# name of dataset
    nrep="numeric",		# number of replications
    datsim="data.frame"	# data (a dataframe with columns: idsim (subject id), irsim (replication index), xsim (simulated x), ysim (simulated response), ydsim (decorrelated replicated data))
  ),
  validity=function(object){
#    cat ("--- Checking NpdeData object ---\n")
    if (length(object@name.simdata)==0) {
      stop ("[ NpdeData : validation ] Please provide a name for the data (dataset or datafile on disk).")
    }
    return(TRUE)
  }
)

###############################
# ECO validity ne semble pas etre appele automatiquement quand on cree un objet => il faut l'appeler dans initialize

setMethod(
  f="initialize",
  signature="NpdeData",
  definition= function (.Object,name.data,name.group,name.predictor, name.response,name.covariates,name.cens,name.miss,name.ipred,units,data){
#    cat ("--- initialising NpdeData Object --- \n")
    if(missing(name.data) || length(name.data)==0) stop ("Please provide a name for the data (dataset or datafile on disk).")
    .Object@name.data<-name.data
    if(missing(name.group)) name.group<-character()
# ECO TODO: reconnaissance automatique (avant affectation a la valeur 2) ?
    if(missing(name.predictor)) name.predictor<-character()
    if(missing(name.response)) name.response<-character()
    if(missing(name.covariates) || length(name.covariates)==0 || name.covariates[1]=="") name.covariates<-character()
    if(missing(name.cens) || length(name.cens)==0 || name.cens=="") name.cens<-character()
    if(missing(name.miss) || length(name.miss)==0 ||name.miss=="") name.miss<-character()
    if(missing(name.ipred) || length(name.ipred)==0 ||name.ipred=="") name.ipred<-character()
    .Object@name.group<-name.group
    .Object@name.predictor<-name.predictor
    .Object@name.response<-name.response
    .Object@name.covariates<-name.covariates
    .Object@name.cens<-name.cens
    .Object@name.miss<-name.miss
    .Object@name.ipred<-name.ipred
    if(missing(units)) units<-list(x="-",y="-")
    if(is.null(units$x)) units$x<-"-"
    if(is.null(units$y)) units$y<-"-"
    ncov<-length(name.covariates)
    if(ncov>0) {
      nunit<-length(units$covariates)
      if(nunit==0) units$covariates<-rep("-",ncov)
      if(nunit>ncov) units$covariates<-units$covariates[1:ncov]
      if(nunit<ncov) {
        length(units$covariates)<-ncov
        units$covariates[(nunit+1):ncov]<-"-"
      }
    }
    .Object@units<-units
    .Object@N<-0
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

setMethod(
  f="initialize",
  signature="NpdeSimData",
  definition= function (.Object,name.simdata){
#    cat ("--- initialising NpdeSimData Object --- \n")
    if(missing(name.simdata)) stop ("Please provide a name for the simulated dataset (dataset or datafile on disk) as a character string.")
    .Object@name.simdata<-name.simdata
    return (.Object )
  }
)

##################################################################################

##' Get/set methods for NpdeData object
##' 
##' Access slots of a NpdeData using the object["slot"] format
##' 
##' @keywords methods
##' @exportMethod [
##' @exportMethod [<-

#### NpdeData
# Getteur
setMethod(
  f ="[",
  signature = "NpdeData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.data"={return(x@name.data)},
    "header"={return(x@header)},
    "sep"={return(x@sep)},
    "na"={return(x@na)},
    "name.group"={return(x@name.group)},
    "name.predictor"={return(x@name.predictor)},
    "name.response"={return(x@name.response)},
    "name.cens"={return(x@name.cens)},
    "name.miss"={return(x@name.miss)},
		"name.covariates"={return(x@name.covariates)},
    "name.ipred"={return(x@name.ipred)},
    "units"={return(x@units)},
    "loq"={return(x@loq)},
    "data"={return(x@data)},
  	"ind"={return(x@ind)},
  				"icens"={return(x@icens)},
  				"not.miss"={return(x@not.miss)},
  	"N"={return(x@N)},
  	"ntot.obs"={return(x@ntot.obs)},
    "nind.obs"={return(x@nind.obs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "NpdeData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.data"={x@name.data<-value},
    "name.group"={x@name.group<-value},
    "name.predictor"={x@name.predictor<-value},
    "name.response"={x@name.response<-value},
    "name.cens"={x@name.cens<-value},
    "name.miss"={x@name.miss<-value},
  	"name.covariates"={x@name.covariates<-value},
  	"name.ipred"={x@name.ipred<-value},
    "units"={x@units<-value},
    "loq"={x@loq<-value},
    "data"={x@data<-value},
  	"ind"={x@ind<-value},	
  				"icens"={x@icens<-value},	
  				"not.miss"={x@not.miss<-value},		
  	"N"={x@N<-value},
  	"ntot.obs"={x@ntot.obs<-value},
    "nind.obs"={x@nind.obs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

#### Simulated data
# Getteur
setMethod(
  f ="[",
  signature = "NpdeSimData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.simdata"={return(x@name.simdata)},
    "nrep"={return(x@nrep)},
    "datsim"={return(x@datsim)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "NpdeSimData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.simdata"={x@name.simdata<-value},
    "nrep"={x@nrep<-value},
    "datsim"={x@datsim<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

##################################################################################
# print/show/showall
# alias in class documentation

##' @S3method print NpdeData
print.NpdeData <- function(x,nlines=10,...) {
    digits<-2;nsmall<-2
    cat("Object of class NpdeData\n")
    cat("    longitudinal data\n")
    if(length(x@name.data)>0)
    cat("Dataset",x@name.data,"\n")
    if(length(x@name.group)>0) {
    st1<-paste(x@name.response," ~ ",paste(x@name.predictor,collapse=" + ")," | ", x@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    cat("    predictor:",x@name.predictor,paste("(",x@units$x,")",sep=""),"\n")
    if(length(x@name.covariates)>0) {
      cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
    }
    if(dim(x@data)[1]>0) {
      if(nlines==0) return()
      cat("Dataset characteristics:\n")
        cat("    number of subjects:    ",x@N,"\n")
        cat("    number of non-missing observations:",x@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(x@nind.obs),digits=digits, nsmall=nsmall), " / ", min(x@nind.obs)," / ",max(x@nind.obs),"\n")
        if(length(x@loq)>0) cat("      LOQ:    ",x@loq,"\n")
#    if(length(x@tab)>0) print(x@tab)
      if(nlines==(-1)) {
        cat("Data:\n")
        print(x@data)
      } else {
        cat("First",nlines,"lines of data:\n")
        nrowShow <- min (nlines , nrow(x@data))
        print(x@data[1:nrowShow,])
      }
    } else cat("No data.\n")
  } else cat("Empty object\n")
}

##' @S3method print NpdeSimData

print.NpdeSimData<-function(x,nlines=10,...) {
    cat("Object of class NpdeSimData\n")
    cat("    dataset:               ",x@name.simdata,"\n")
    if(length(x@nrep)>0) cat("    number of replications:",x@nrep,"\n") else cat("    x currently empty\n")
    if(length(x@datsim)>0) {
      if(nlines==0) return()
      if(nlines==(-1)) {
        cat("Data:\n")
        print(x@datsim)
      } else {
        cat("First",nlines,"lines of simulated data:\n")
        nrowShow <- min (nlines , nrow(x@datsim ))
        print(x@datsim[1:nrowShow,])
      }
    } else cat("No data.\n")
}

setMethod("show","NpdeData",
# show.NpdeData<-
	function(object) {
    cat("Object of class NpdeData\n")
    if(length(object@ntot.obs)==0) cat("    no data\n") else {
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictor,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(object@name.covariates)>0) cat("    Covariates:",object@name.covariates,"\n")
      cat("This object has the following components:\n") 
      cat("     data: data\n")
      cat("     with",object@N,"subjects\n")
      cat("     ",object@ntot.obs,"observations\n")
      cat("The data has the following components\n")
      cat("     X:",object@name.predictor,"\n")
      cat("     Y:",object@name.response,"\n")
      if(length(object@name.ipred)>0) cat("     individual model predictions:", object@name.ipred,"\n")
      if(length(object@name.miss)>0) cat("     missing data:",object@name.miss," (1=missing)\n")
      if(length(object@name.cens)>0) cat("     censored data:",object@name.cens," (1=censored)\n")
      if(length(object@loq)>0) cat("      LOQ:    ",object@loq,"\n")
#      cat("     \n")
  }
  }
)

#show.NpdeSimData<-
setMethod("show","NpdeSimData",
	function(object) {
    cat("Object of class NpdeSimData\n")
    cat("    dataset:               ",object@name.simdata,"\n")
    if(length(object@nrep)>0) cat("    number of replications:",object@nrep,"\n") else cat("    object currently empty\n")
  }
)


#
#' Brief summary of an object
#'
#' Prints a brief summary of an object
#' 
#' @name showall
#' @aliases showall showall.NpdeData showall,NpdeData-method showall.NpdeSimData showall,NpdeSimData-method  showall.NpdeRes  showall.NpdeObject
#' @param object an object
#' @keywords print
##' @S3method showall NpdeData
##' @export
# Could be print, with only head of data

#setMethod("show","NpdeSimData",

#setMethod("showall","NpdeData",
showall.NpdeData <- function(object) {
    digits<-2;nsmall<-2
    cat("Object of class NpdeData\n")
    cat("    longitudinal data\n")
    if(length(object@name.data)>0) {
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictor,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    cat("    subject identifier:    ",object@name.group,"\n")
    cat("    predictor:       ",object@name.predictor, paste("(",object@units$x,")",sep=""),"\n")
    cat("    response:        ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
    if(length(object@name.covariates)>0) {
    	cat("    covariates:",paste(paste(object@name.covariates," (", object@units$covariates,")",sep=""),collapse=", "),"\n")
    }
    cat("This object has the following components:\n") 
    cat("     data: data\n")
    cat("     with",object@N,"subjects\n")
    cat("     ",object@ntot.obs,"observations\n")
    cat("The data has the following components\n")
    cat("     X:",object@name.predictor,"\n")
    cat("     Y:",object@name.response,"\n")
    if(length(object@name.ipred)>0) cat("     individual model predictions:", object@name.ipred,"\n")
    if(length(object@name.miss)>0) cat("     missing data:",object@name.miss," (1=missing)\n")
    if(length(object@name.cens)>0) cat("     censored data:",object@name.cens," (1=censored)\n")
    if(length(object@loq)>0) cat("      LOQ:    ",object@loq,"\n")
    cat("Dataset characteristics:\n")
    cat("    number of subjects:    ",object@N,"\n")
    if(object@N>0) {
      cat("    number of non-missing observations:",object@ntot.obs,"\n")
      cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@orig)>0) print(object@orig)
    }
    if(dim(object@data)[1]>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@data))
      print(object@data[1:nrowShow,])
    } else cat("No data.\n")
  } else cat("Empty object\n")
    }
#)

##' @S3method showall NpdeSimData

#setMethod("showall","NpdeSimData",
showall.NpdeSimData<-
					function(object) {
						cat("Object of class NpdeSimData\n")
						cat("    dataset:               ",object@name.simdata,"\n")
						if(length(object@nrep)>0) cat("    number of replications:",object@nrep,"\n") else cat("    object currently empty\n")
					}
#)

##' @S3method summary NpdeData

summary.NpdeData <- function(object, print=TRUE, ...) {
	  if(length(object@data)==0) {
			cat("Object of class NpdeData, empty.\n")
			return()
		}
	  res<-list(N=object@N,data=object@data, ntot.obs=object@ntot.obs,nind.obs=object@nind.obs)
	  if(length(object@loq)>0) res$loq<-object@loq
	  invisible(res)
}

##' @S3method subset NpdeData

subset.NpdeData<-function (x, subset, ...) {
    if (missing(subset)) 
        return(x)
    else {
        e <- substitute(subset)
	xdat<-x["data"]
        r <- eval(e, xdat, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    x1<-x
    x1["data"]<-x["data"][r,,drop=FALSE]
    if(length(x1["not.miss"])>0) {
      x1["not.miss"]<-x["not.miss"][r]
      x1["icens"]<-which(!x1["not.miss"])
    }
    id<-x1["data"][,x1["name.group"]]
    x1["N"]<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-c(nind.obs[match(unique(id),names(nind.obs))])
    x1["nind.obs"]<-nind.obs
    x1["ntot.obs"]<-length(id)
    x1["ind"]<-rep(1:x1["N"],times=nind.obs)
    return(x1)
}

##################################################################################

##' @S3method set.plotoptions NpdeData

set.plotoptions.NpdeData<-function(object, ...) {
# setting default plot options
  plot.opt<-list(
# Options for plot types
    ilist=c(1:object["N"]),
    level=0:1,
    plot.loq=TRUE,			# Whether data under the LOQ should be plotted
    line.loq=TRUE,			# Whether an horizontal line should be plotted at Y=LOQ
    # General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main="",				# title
    sub="",
    xlab="",
    ylab="",
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1,
    col="black",
    lcol="black",
    col.pobs="black",		# color to plot observed/censored data
    col.pcens="red",
    col.lobs="black",		# color to plot observed/censored data
    pch=20,
    pch.pobs=20,				# symbol for observed/censored data
    pch.pcens=8,
    ablinecol="black",
    ablinelty=2,
    ablinelwd=2
    )
        
     if(is.null(plot.opt$name.X))
        plot.opt$name.X<-object["name.predictor"][1]
    plot.opt$xlab<-paste(plot.opt$name.X," (",plot.opt$units$x,")", sep="")
     if(length(object["name.response"])>0)
    plot.opt$ylab<-paste(object["name.response"]," (",object["units"]$y,")", sep="")
   return(plot.opt)
}



#' Replace graph options
#'
#' This function is used to replace graph options (available in the prefs slot of the NpdeData object) for plotting NpdeData objects
#' 
#' @usage replace.plotoptions(plot.opt,...)
#' @param plot.opt a list of graphical preferences
#' @param \dots names and values of the options to be replaced
#' 
#' @return an updated list of options for graphs
#' @details See documentation for a list of available options. During replacement, invalid (including misspelled) options will raise warnings. 
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @keywords plot internal

replace.plotoptions<-function(plot.opt,...) {
	args1<-match.call(expand.dots=TRUE)
# These arguments are used by other functions and may be passed on via "...", so we want to ignore them. Other arguments not in list will raise warnings
	legacy<-c("plot.type","namsav","namgr","loq")
	if(length(args1)>2) {
		# General arguments: col, pch
		i1<-match("col",names(args1))
		if(!is.na(i1)) {
			plot.opt$col<-eval(args1[[i1]])
			if(is.na(match("pcol",names(args1)))) plot.opt$pcol<-eval(args1[[i1]])
			if(is.na(match("lcol",names(args1)))) plot.opt$lcol<-eval(args1[[i1]])
			if(is.na(match("ablinecol",names(args1)))) plot.opt$ablinecol<-eval(args1[[i1]])
			if(is.na(match("col.pobs",names(args1)))) plot.opt$col.pobs<-eval(args1[[i1]])
			if(is.na(match("col.lobs",names(args1)))) plot.opt$col.lobs<-eval(args1[[i1]])
			if(is.na(match("col.pcens",names(args1)))) plot.opt$col.pcens<-eval(args1[[i1]])
#			if(is.na(match("col.lcdf",names(args1)))) plot.opt$col.lcdf<-eval(args1[[i1]])
#			if(is.na(match("col.fill",names(args1)))) plot.opt$col.fill<-eval(args1[[i1]])
#			if(is.na(match("col.fillcdf",names(args1)))) plot.opt$col.fillcdf<-eval(args1[[i1]])
				}
		i1<-match("lty",names(args1))
		if(!is.na(i1)) {
			plot.opt$lty<-eval(args1[[i1]])
			if(is.na(match("lty.lobs",names(args1)))) plot.opt$lty.lobs<-eval(args1[[i1]])
			if(is.na(match("ablinelty",names(args1)))) plot.opt$ablinelty<-eval(args1[[i1]])
		}
		i1<-match("lwd",names(args1))
		if(!is.na(i1)) {
			plot.opt$lwd<-eval(args1[[i1]])
			if(is.na(match("lwd.lobs",names(args1)))) plot.opt$lwd.lobs<-eval(args1[[i1]])
			if(is.na(match("ablinelwd",names(args1)))) plot.opt$ablinelwd<-eval(args1[[i1]])
		}
		i1<-match("pch",names(args1))
		if(!is.na(i1)) {
			plot.opt$pch<-eval(args1[[i1]])
			if(is.na(match("pch.pobs",names(args1)))) plot.opt$pch.pobs<-eval(args1[[i1]])
			if(is.na(match("pch.pcens",names(args1)))) plot.opt$pch.pcens<-eval(args1[[i1]])
		}
		# Other arguments
		for(i in 3:length(args1)) {
			if(match(names(args1)[i],names(plot.opt),nomatch=0)>0) {
				#    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
				if(!is.null(eval(args1[[i]]))) plot.opt[[names(args1)[i]]]<-eval(args1[[i]])
				} else {
					if(is.na(match(names(args1)[i],legacy))) cat("Argument",names(args1)[i],"not available, check spelling.\n")
				}
			}
		}
	return(plot.opt)
	}

#' Plots a NpdeData object
#'
#' Plots the data in a NpdeData object
#' 
#' @param x a NpdeData object
#' @param y unused, here for compatibility with the base plot function
#' @param \dots additional graphical parameters to be passed on to the plot
#' @details The default plot is a spaghetti plot of all the data, with a line joining the observations for each subject. If censored data is present, it is shown with a different symbol and colour.
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
##' @seealso \code{\link{set.plotoptions}}
#' @keywords plot
#' @examples
#' 
#' data(theopp)
#' 
#' x<-npdeData(theopp,name.group="ID",name.predictor="Time",name.response="Conc", 
#' name.covariates=c("Wt"),units=list(x="hr",y="mg/L",covariates="kg"))
#' plot(x)
#' 
##' @importFrom graphics plot
##' @method plot NpdeData
##' @export

plot.NpdeData <- function(x, y, ...) {
# Plot the data, either as points or as lines grouped by x@name.group
    if(length(x@data)>0) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"b"
    plot.opt<-set.plotoptions(x)
    plot.opt$new<-TRUE
    plot.opt$xlab<-paste(x@name.predictor," (",x@units$x,")",sep="")
    plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
    plot.opt<-replace.plotoptions(plot.opt,...)
    logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
    if(plot.opt$new) par(mfrow=c(1,1))
    tab<-x@data[x@ind,] # remove missing data
    has.cens<-(length(x@name.cens)>0)
    if(has.cens && max(tab[,x@name.cens])==0) has.cens<-FALSE
      if(plot.type=="p" | plot.type=="b") {
      	if(has.cens) {
        plot(tab[,x@name.predictor],tab[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,sub=plot.opt$sub,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab,type="n")
        points(tab[tab[,x@name.cens]==0,x@name.predictor], tab[tab[,x@name.cens]==0,x@name.response],col=plot.opt$col.pobs, pch=plot.opt$pch.pobs,cex=plot.opt$cex)
        points(tab[tab[,x@name.cens]==1,x@name.predictor], tab[tab[,x@name.cens]==1,x@name.response],col=plot.opt$col.pcens, pch=plot.opt$pch.pcens,cex=plot.opt$cex)
      	} else 
      		plot(tab[,x@name.predictor],tab[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col.pobs,pch=plot.opt$pch.pobs,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,sub=plot.opt$sub,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
      }
      if(plot.type=="l") {
        plot(tab[,x@name.predictor],tab[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main,sub=plot.opt$sub, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
      }
      if(plot.type=="l" | plot.type=="b") {
        for(isuj in unique(tab[,x@name.group])) {
          lines(tab[tab[,x@name.group]==isuj,x@name.predictor], tab[tab[,x@name.group]==isuj,x@name.response],col=plot.opt$col.lobs, lty=plot.opt$lty,lwd=plot.opt$lwd)
      }
      }
    if(has.cens) abline(h=x@loq,col=plot.opt$ablinecol,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)

    } else cat("No data to plot.\n")
  }

##################################################################################
# @usage read.npdeData(object, header, sep, na.strings, detect)

#' Read observed data
#' 
#' Creates an NpdeData object containing the observed data, either from disk or from a dataframe
#' 
#'
#' @param header boolean indicating whether the file has a header (mandatory if 
#' detect is TRUE)
#' @param sep field separator (for files on disk)
#' @param na.strings strings to be considered as indicating NA
#' @param detect a boolean; if TRUE, automatic recognition of names will be attempted to detect necessary items (longitudinal data structure, missing data and censoring information)
#' @return an object of class \code{"\linkS4class{NpdeData}"}
#' @keywords internal

setMethod("read.npdeData","NpdeData",
  function(object,header,sep,na.strings,detect,verbose) {
    ow <- options("warn")
    options("warn"=-1)
# ce test devrait aller dans la definition de la classe
    if(class(object@name.data)!="character") {
    cat("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
    return("Creation of npdeData failed")
  }
    if(exists(object@name.data)) {
      if(verbose) cat("Using the object called",object@name.data,"in this R session as the data.\n")
      dat<-get(object@name.data)
    } else {
      if(verbose) cat("Reading data from file",object@name.data,"\n")
      dat<-try(read.table(object@name.data,header=header,sep=sep,na.strings=na.strings))
      if(class(dat)=="try-error") stop("The file ",object@name.data," does not exist. Please check the name and path.\n")      
      if(verbose) {
      	cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
      	print(head(dat))
      }
    }
    if(dim(dat)[2]<2) {
      cat("The dataset contains only one column. To compute npde, we need at least 3 columns, with subject ID, predictor (at least one) and response. \nPlease check the field separator, currently given as:", paste("sep=\"",sep,"\"",sep=""), "\n")
      return("Creation of npdeData failed")
    }
# Automatic recognition of columns 
#    ID (one of id, subject or sujet regardless of case)
#    response (one of Y, conc, concentration, resp, response regardless of case)
#    predictors (time and/or dose, regardless of case)
# ECO TODO: improve automatic recognition ?
# check that we have at least a column id, response and X
    if(!is.na(as.integer(object@name.group))) {
# group given as a column number
      object@name.group<-colnames(dat)[as.integer(object@name.group)]
    }
    if(is.na(object@name.group) || object@name.group=="") {
    	if(!detect) {
    		cat("Missing ID column and automatic detection is OFF. Please provide a valid name for the ID column\n")
    		return("Creation of npdeData failed")
    	}
      if(verbose) cat("Missing ID column, attempting to detect it\n")
      object@name.group<-""
      i1<-match("id",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) {
      	i1<-c(match(c("subject","sujet","group"),tolower(colnames(dat))))
      }
      if(length(i1)>0) {
        object@name.group<-colnames(dat)[i1[1]]
        if(verbose) cat("    no name for the group variable (ID) given, will use column --",object@name.group,"-- in the dataset.\n")
      }
    }
    if(object@name.group=="" | is.na(match(object@name.group,colnames(dat)))) {
      cat("Please provide a name for the ID column.\n")
      return("Creation of npdeData failed")
    }
# Predictors
   i1<-as.integer(object@name.predictor[!is.na(as.integer(object@name.predictor))])
    if(length(i1)>0) { 
      object@name.predictor[!is.na(as.integer(object@name.predictor))]<- colnames(dat)[i1]
    }
    if(is.na(object@name.predictor) | length(object@name.predictor)==0 | (length(object@name.predictor)==1 & object@name.predictor[1]=="")) {
    	if(!detect) {
    		cat("Missing X column and automatic detection is OFF. Please provide a valid name for the column with the predictor.\n")
    		return("Creation of npdeData failed")
    	}
    	if(verbose) cat("Missing predictor column, attempting to detect it\n")
      object@name.predictor<-""
      i1<-c(match(c("xobs","time","temps","tps","tim","x","dose"), tolower(colnames(dat))))
      i1<-i1[!is.na(i1)]
      if(length(i1)>0) {
        object@name.predictor<-colnames(dat)[i1][1]
        if(verbose) cat("    no name for the predictor variable given, will use column(s) --",object@name.predictor,"-- in the dataset.\n")
      }
    }
    id1<-match(object@name.predictor,colnames(dat),nomatch=0)
    if(length(id1[id1==0])>0) {
    	if(verbose) cat("    cannot find column(s) --",object@name.predictor[id1==0],"-- dropping them from the data.\n")
    }
    xnam<-object@name.predictor[id1>0]
    if(length(xnam)==0) object@name.predictor<-"" else object@name.predictor<-xnam
    if(length(xnam)==0) {
      cat("Please provide at least one predictor.\n")
      return("Creation of npdeData failed")
    }
# Response
    if(!is.na(as.integer(object@name.response))) { 
# response given as a column number
      object@name.response<-colnames(dat)[as.integer(object@name.response)]
    }
    if(is.na(object@name.response) || object@name.response=="") {
    	if(!detect) {
    		cat("Missing response column and automatic detection is OFF. Please provide a valid name for the column with the response.\n")
    		return("Creation of npdeData failed")
    	}
    	if(verbose) cat("Missing response column, attempting to detect it\n")
      object@name.response<-""
      i1<-match("y",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) { 
        i1<-c( match(c("yobs","resp","conc"),tolower(colnames(dat))), grep("response",tolower(colnames(dat)),fixed=TRUE),grep("concentration", tolower(colnames(dat)),fixed=TRUE))
        i1<-i1[!is.na(i1)]
      }
      if(length(i1)>0) {
        object@name.response<-colnames(dat)[i1[1]]
        if(verbose) cat("    no name for the response variable given, will use column --",object@name.response,"-- in the dataset.\n")
      }
    }
    if(is.na(object@name.response)) object@name.response<-""
    if(object@name.response=="" | is.na(match(object@name.response,colnames(dat)))) {
      cat("Please provide a name for the response column.\n")
      return("Creation of npdeData failed")
    }
# ECO TODO: verifier que les colonnes existent et sinon corriger
    
# IPRED : column with individual predictions
    detect.ipred<-FALSE
    if(length(object@name.ipred)>0 && !is.na(as.integer(object@name.ipred))) # ipred given as a column number
    	object@name.ipred<-colnames(dat)[as.integer(object@name.ipred)]
    if(length(object@name.ipred)>0 && match(object@name.ipred,colnames(dat),nomatch=0)==0) {
    	if(detect & verbose) cat("Can't find a column named",object@name.ipred,"in the dataset for individual predictions, will attempt automatic detection.\n")
    	object@name.ipred<-character()
    }
    if(length(object@name.ipred)==0 || is.na(object@name.ipred)) detect.ipred<-TRUE
    if(detect.ipred) {
    	i1<-c(grep("ipred",tolower(colnames(dat)),fixed=T))
    	if(length(i1)>0) {
    		object@name.ipred<-colnames(dat)[i1[1]]
    		if(detect.ipred & verbose) cat("    assuming that individual predictions are given in column --",object@name.ipred,"-- in the dataset (to ignore this column, add the argument detect=FALSE in the call to npdeData()).\n")
    	}
    }
# CENS : column indicating censoring
    detect.cens<-FALSE
    if(length(object@name.cens)>0 && !is.na(as.integer(object@name.cens))) # cens given as a column number
    	object@name.cens<-colnames(dat)[as.integer(object@name.cens)]
    if(length(object@name.cens)>0 && match(object@name.cens,colnames(dat),nomatch=0)==0) {
    	if(detect & verbose) cat("Can't find a column named",object@name.cens,"in the dataset containing censoring, will attempt automatic detection.\n")
    	object@name.cens<-character()
    }
    if(length(object@name.cens)==0 || is.na(object@name.cens)) detect.cens<-TRUE
    if(detect.cens) {
    	i1<-c(grep("cens",tolower(colnames(dat)),fixed=T))
    	if(length(i1)>0) {
    		object@name.cens<-colnames(dat)[i1[1]]
    		if(detect.cens & verbose) cat("    assuming that censoring information is given in column --",object@name.cens,"-- in the dataset (to ignore this column, add the argument detect=FALSE in the call to npdeData()).\n")
    	}
    }
    if(length(object@name.cens)>0) { # checking validity of censoring column
    	if(!isTRUE(all.equal(sort(unique(dat[,object@name.cens]), na.last=TRUE),as.integer(c(0,1))))) {
    		if(verbose) cat("The column with censoring information should only contain 0 and 1s.\n")
    	  object@name.cens<-character()
    }}
# MDV : column indicating missing data
    detect.miss<-FALSE
    if(length(object@name.miss)>0 && !is.na(as.integer(object@name.miss))) # miss given as a column number
    	object@name.miss<-colnames(dat)[as.integer(object@name.miss)]
    if(length(object@name.miss)>0 && match(object@name.miss,colnames(dat),nomatch=0)==0) {
    	if(detect & verbose) cat("Can't find a column named",object@name.miss,"in the dataset containing missing data status, will attempt automatic detection.\n")
    	object@name.miss<-character()
    }
    if(length(object@name.miss)==0 || is.na(object@name.miss)) detect.miss<-TRUE
    if(detect.miss) {
    	i1<-c(grep("mdv",tolower(colnames(dat)),fixed=T), grep("miss",tolower(colnames(dat)),fixed=T))
    	if(length(i1)>0) {
    		object@name.miss<-colnames(dat)[i1[1]]
    		if(detect.miss & verbose) cat("    assuming that column --",object@name.miss,"-- in the dataset contains missing data information (to ignore this column, add the argument detect=FALSE in the call to npdeData()).\n")
    	}
    }
    if(length(object@name.miss)>0) { # checking validity of MDV column
    	if(!isTRUE(all.equal(sort(unique(dat[,object@name.miss]), na.last=TRUE),as.integer(c(0,1)))) & !isTRUE(all.equal(sort(unique(dat[,object@name.miss]), na.last=TRUE),as.integer(c(0)))) & !isTRUE(all.equal(sort(unique(dat[,object@name.miss]), na.last=TRUE),as.integer(c(1))))) {
    		if(verbose) cat("The column with information about missing data should only contain 0 and 1s.\n")
    		object@name.miss<-character()
    	}}
# Covariates
    if(length(object@name.covariates)>0 & object@name.covariates[1]!="") {
      i1<- as.integer(object@name.covariates[!is.na(as.integer(object@name.covariates))])
      object@name.covariates[!is.na(as.integer(object@name.covariates))]<- colnames(dat)[i1]
    }
# Saving covariates in the original format in ocov, transforming binary covariates in dat to factors - No: here we can keep covariates as factors
#     object@ocov<-dat[,object@name.covariates,drop=FALSE]
#     for(icov in object@name.covariates) {
#     	if(length(unique(dat[,icov]))==2) dat[,icov]<-as.integer(factor(dat[,icov]))-1
#     }   
    
    if(nchar(object@name.group)*length(object@name.predictor)* nchar(object@name.response)<=0) {
      stop("Please check the structure of the data file and provide information concerning which columns specify the group structure (ID), the predictors (eg dose, time) and the response (eg Y, conc). See documentation for automatic recognition of column names for these elements.\n")
    }
# Data
    all.names<-c(object@name.group,object@name.predictor,object@name.response, object@name.covariates,object@name.miss,object@name.cens,object@name.ipred)
    tab<-dat[,all.names]
# Index (ID may be numbers, strings,...)
    id<-tab[,object@name.group]
    # ECO TODO: et si un sujet n'a que des donnees NA ???
    object@N<-length(unique(id))
    nind.obs.full<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs.full<-nind.obs.full[match(unique(id),names(nind.obs.full))]
    tab<-data.frame(index=rep(1:object@N,times=nind.obs.full),tab)
# Missing data
    if(length(object@name.miss)>0) mdv<-tab[,object@name.miss] else {
    	mdv<-rep(0,length(id))
    	object@name.miss<-"mdv"
    } 
    mdv[is.na(tab[,object@name.response])]<-1
    tab[,object@name.miss]<-mdv    
    object@data<-tab
    object@ind<-which(mdv==0)
    icens<-numeric()
    if(length(object@name.cens)>0) {
    	icens<-which(mdv==0 & dat[,object@name.cens]==1)
    } 
    object@icens<-icens
    object@not.miss<-(mdv==0)
# ECO TODO: what about missing data in covariates & predictor columns
    if(length(object@name.covariates)>0 && sum(is.na(object@data[object@not.miss,object@name.covariates]))>0) {
    	tab<-object@data
    	for(icov in object@name.covariates) {
    		for(i in 2:dim(tab)) {
    			if(is.na(tab[i,icov])) tab[i,icov]<-tab[(i-1),icov]    			
    		}    		
    	}
    	object@data<-tab
    }
#    for(i in object@name.covariates) 
#      object@data[is.na(object@data[,i]),object@name.miss]<-1
    tb1<-tab[tab[,object@name.miss]==0,]
    id1<-tb1[,1]
    object@ntot.obs<-dim(tb1)[1] # total number of observations
    nind.obs<-tapply(id1,id1,length) # individual numbers of observations (1xN)
    nind.obs<-nind.obs[match(unique(id1),names(nind.obs))]
    object@nind.obs<-c(nind.obs)    

#    object@names<-list(group=object@name.group,predictors=object@name.predictor, response=object@name.response, covariates=object@name.covariates)
    options(ow) # reset
    validObject(object)
    return(object)
  }
)

#' Read simulated data
#' 
#' Creates an NpdeSimData object containing the simulated data, either from disk or from a dataframe
#' 
#' @param header boolean indicating whether the file has a header (mandatory if 
#' detect is TRUE)
#' @param sep field separator (for files on disk)
#' @param na.strings strings to be considered as indicating NA
#' @return an object of class NpdeSimData
#' @keywords internal

setMethod("read.npdeSimData","NpdeSimData",
	function(object,header,sep,na.strings,verbose) {
    ow <- options("warn")
    options("warn"=-1)
    if(exists(object@name.simdata)) {
    	if(verbose) cat("Using the object called",object@name.simdata,"in this R session as the data.\n")
      dat<-get(object@name.simdata)
    } else {
    	if(verbose) cat("Reading data from file",object@name.simdata,"\n")
      if(missing(header)) header<-TRUE
      if(missing(sep)) sep<-""
      if(missing(na.strings)) na.strings<-c(".","NA")
      dat<-try(read.table(object@name.simdata,header=header,sep=sep,na.strings=na.strings))
      if(class(dat)=="try-error") stop("The file ",object@name.simdata," does not exist. Please check the name and path.\n")      

      x1<-unlist(strsplit(as.character(dat[1,1]),",",fixed=TRUE))
      mysep<-""
      if(length(x1)>1) mysep<-","
      x1<-unlist(strsplit(as.character(dat[1,1]),";",fixed=TRUE))
      if(length(x1)>1) mysep<-";"
      if(mysep!="") dat<-read.table(object@name.simdata,na.strings=c(".","NA"),sep=mysep)
      if(!is.numeric(dat[1,1]))
        dat<-read.table(object@name.simdata,na.strings=c(".","NA"),header=TRUE,sep=mysep)
      if(!is.numeric(dat[1,1])) {
        cat("The format of the file containing the simulated data is unknown.\n")
        cat("Please use a standard R table format, with or without header,\n")
        cat("and with one of the following separators: \n")
        cat("         TAB or space(s), commas (',') or semicolons (';')\n")
        cat("Also note that a dot should be used to indicate digits in numbers.\n")
        stop("Exiting npde\n")
      }
    	if(verbose) {
    		cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
        print(head(dat))
    	}
    }
    colnames(dat)<-c("idsim","xsim","ysim")
    object@datsim<-dat

    options(ow) # reset
    validObject(object)
    return(object)
  }
)

##################################################################################
#
#' Creates a NpdeData object
#'
#' This function is used to create a NpdeData object, representing a longitudinal data structure, and fill it with data from a dataframe or a file on disk
#' 
#' @usage npdeData(name.data,header=TRUE,sep="",na.strings=c(".","NA"),name.group, name.predictor,name.response, name.covariates,name.cens,name.miss,name.ipred, units=list(x="",y="",covariates=c()),detect=TRUE,verbose=FALSE)
#' @param name.data name of the file containing the observed data, or a dataframe
#' containing the observed data
#' @param sep field separator (for files on disk)
#' @param na.strings strings to be considered as indicating NA
#' @param header boolean indicating whether the file has a header (mandatory if 
#' detect is TRUE)
#' @param name.group name/number of the column in the observed data containing the 
#' patient ID (if missing and detect is TRUE, columns named id, subject or sujet 
#' (regardless of case) will be assumed to contain this information)
#' @param name.predictor name/number of the column in the observed data containing 
#' the independent variable X (if missing and detect is TRUE, columns named xobs, 
#' time, dose, x, temps, tim (regardless of case) will be assumed to 
#' contain this information)
#' @param name.response name/number of the column in the observed data containing 
#' the dependent variable Y (if missing and detect is TRUE, columns named yobs, 
#' response, resp, conc, concentration (regardless of case) will be assumed to 
#' contain this information)
#' @param name.miss name/number of the column containing information about missing 
#' data (MDV) (if missing and detect is TRUE, column called mdv or miss 
#' (regardless of case) will be assumed to contain this information)
#' @param name.cens name/number of the column containing information about censored 
#' data (cens) (if missing and detect is TRUE, column with a name containing cens 
#' (regardless of case) will be assumed to contain this information)
#' @param name.covariates name/number of the column(s) containing covariate 
#' information (optional)
#' @param name.ipred name/number of the column(s) with individual predictions
#' (ipred)  (if missing and detect is TRUE, column with a name containing ipred 
#' (regardless of case) will be assumed to contain this information)
#' @param units a list with components x, y and cov (optional), specifying the
#' units respectively for the predictor (x), the response (y), and the covariates 
#' (a vector of length equal to the number of covariates). Units will default to (-) if not given.	
#' @param detect a boolean controlling whether automatic recognition of columns in the dataset is on, defaults to TRUE
#' @param verbose whether to print warning messages, defaults to FALSE (set to TRUE to check how data is being handled)
#' 
#' @return an object of class NpdeData
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords models
#' @export
#' @examples
#' 
#' data(theopp)
#' 
#' x<-npdeData(theopp) # Automatic detection
#' print(x)
#' x<-npdeData(theopp,name.group="ID",name.predictor="Time",name.response="Conc", 
#' name.covariates=c("Wt"),units=list(x="hr",y="mg/L",covariates="kg")) # Explicit
#' print(x)
#' plot(x)

npdeData<-function(name.data,header=TRUE,sep="",na.strings=c(".","NA"),name.group, name.predictor,name.response, name.covariates,name.cens,name.miss,name.ipred, units=list(x="",y="",covariates=c()),detect=TRUE,verbose=FALSE) {
# setting proper types for the NpdeData class
  if(missing(name.data) ||length(name.data)==0) {
    cat("Error in npdeData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of NpdeData failed")
  }
  if(is.data.frame(name.data)) name.data<-deparse(substitute(name.data))
  if(missing(name.group)) name.group<-"" else name.group<-as.character(name.group)
  if(missing(name.predictor)) name.predictor<-"" else name.predictor<-as.character(name.predictor)
  if(missing(name.response)) name.response<-"" else  name.response<-as.character(name.response)
  if(missing(name.covariates) || name.covariates[1]==0) name.covariates<-character() else name.covariates<-as.character(name.covariates)
  if(missing(name.miss) || name.miss==0) name.miss<-character() else name.miss<-as.character(name.miss)
  if(missing(name.cens) || name.cens==0) name.cens<-character() else name.cens<-as.character(name.cens)
  if(missing(name.ipred) || name.ipred==0) name.ipred<-character() else name.ipred<-as.character(name.ipred)
  if(missing(detect)) detect<-TRUE
  x<-new(Class="NpdeData",name.data=name.data,name.group=name.group, name.predictor=name.predictor,name.response=name.response, name.covariates=name.covariates,name.cens=name.cens,name.miss=name.miss, name.ipred=name.ipred,units=units)
#  showall(x)
  if(detect & verbose) cat("Automatic detection of variables is ON. The program will attempt to detect both mandatory variables (ID, X, Y) and optional variables (IPRED, MDV, CENS) when they are not specifically given or when the user-specified names are not found in the dataset, by looking in the names of the columns (to override this behaviour, please use argument detect=FALSE in the call to npdeData().\n")
  x1<-read.npdeData(x,header=header,sep=sep,na.strings=na.strings,detect=detect, verbose=verbose)
  if(class(x1)!="character") {
  if(length(x1["name.cens"])==0) loq<-NA else {
  if(sum(x1["data"][x1["data"][,x1["name.miss"]]==0,x1["name.cens"]])>0) {
    yloq<-x1["data"][x1["data"][,x1["name.cens"]]==1 & x1["data"][,x1["name.miss"]]==0,x1["name.response"]]
    if(length(unique(yloq))==1) {
      loq<-unique(yloq)
      if(verbose) cat("Same LOQ for all missing data, loq=",loq,"\n")
    } else {
      loq<-min(unique(yloq),na.rm=TRUE)
      if(verbose)cat("There are different LOQ for different observations, setting loq to the lowest value of",loq,"\n")
    }
    }
    x1["loq"]<-loq
  }
  if(verbose) {
  	cat("\n\nThe following NpdeData object was successfully created:\n\n")
  	print(x1,nlines=0)
  }
  } else x1<-"Creation of NpdeData failed"
  return(x1)
}

npdeSimData<-function(npde.data,name.simdata,verbose=FALSE) {
	if(is.data.frame(name.simdata)) name.simdata<-deparse(substitute(name.simdata))
	ierror<-FALSE
	if(missing(npde.data)) {
		ierror<-TRUE
		cat("   Error: Missing first argument.\n")
	}
	if(!ierror) {
	  x1<-try(class(npde.data))
	  if(class(x1)=="try-error") {
	  	ierror<-TRUE
	  	cat("   Error:", deparse(substitute(npde.data)),"does not exist.\n")
	  }
	  if(!ierror && x1!="NpdeData") {
	  	ierror<-TRUE 
	  	cat("   Error:", deparse(substitute(npde.data)),"is not a NpdeData object.\n")
	  }
	}
	if(ierror) {
			cat("Function npdeSimData requires two mandatory arguments: first, a NpdeData object created by a call to npdeData() (see help page for the syntax of that function), and the name of a matching dataset containing the simulated data (either a file on disk or a data.frame. Please refer to the documentation for details and examples.\n")
		return("Creation of NpdeSimData failed")
	}
	x1<-new(Class="NpdeSimData",name.simdata=name.simdata)
	x<-read.npdeSimData(x1,verbose=verbose)
	if(sum(npde.data["data"][,npde.data["name.miss"]])>0) {
		if(verbose) cat("There are rows with MDV=1 in the original dataset, the corresponding rows will be removed from the simulated dataset.\n")
	}
	nrep<-dim(x@datsim)[1]/dim(npde.data@data)[1]
	x@nrep<-as.integer(nrep)
	if(nrep<1000 & verbose) {
		cat("Warning: the number of simulations is",nrep,"which may be too small.\n")
		cat("We advise performing at least 1000 simulations to compute npde.\n")
	} 
	irsim<-rep(1:nrep,each=dim(npde.data@data)[1])
	x@datsim$irsim<-irsim
	return(x)
}    
##################################################################################
