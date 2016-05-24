##################################################################################

##' Class "NpdeObject"
##' 
##' An object of class NpdeObject
##' 
##' @name NpdeObject-class
##' @aliases NpdeObject NpdeObject-class, show,NpdeObject-method print,NpdeObject-method 
##' showall,NpdeObject-method summary,NpdeObject-method test,NpdeObject-method
##' [,NpdeObject-method [<-,NpdeObject-method
##' npde.main,NpdeObject npde.save,NpdeObject npde.graphs,NpdeObject plot,NpdeObject
##' @docType class
##' @section Objects from the Class: NpdeObject objects are typically created by calls to \code{\link{npde}} or \code{\link{autonpde}}. They contain the following slots:
##' 
##' \describe{
##' \item{data}{an object of class NpdeData, containing the observed data}
##' \item{sim.data}{an object of class NpdeSimData, containing the simulated data}
##' \item{res}{an object of class NpdeRes, containing the results}
##' \item{options}{a list of options}
##' \item{prefs}{a list of graphical preferences for the plots}
##' }
##' @section Methods:
##' \describe{
##'   \item{print(x):}{Prints a summary of object}
##'   \item{show(x):}{Prints a short summary of object}
##'   \item{showall(x):}{Prints a detailed summary of object}
##'   \item{plot(x):}{Diagnostic and other plots. More details can be found in \code{\link{plot.NpdeObject}}}
##'   \item{summary(x):}{Returns a summary of object x in list format}
##'   \item{gof.test(x, which="npde", parametric=TRUE, ...):}{Returns goodness-of-fit tests}
##'   \item{set.plotoptions(x):}{Sets options for graphs (internal method used in plots)}
##' }
##' @seealso \code{\link{npde}}, \code{\link{autonpde}}, \code{\link{NpdeData}}, \code{\link{NpdeSimData}}, \code{\link{NpdeRes}}, \code{\link{gof.test}}
##' @keywords classes
##' @examples
##' 
##' methods(class="NpdeObject")
##' 
##' showClass("NpdeObject")
##' 
##' @exportClass NpdeObject

setClass(Class="NpdeObject",
  representation=representation(
    data="NpdeData",		# Data
    results="NpdeRes",		# Fit results
    sim.data="NpdeSimData", 	# Simulated data
    options="list",		# Options and parameters for algorithm
    prefs="list"		# Options for graphs
  ),
  validity=function(object){
#    cat ("--- Checking NpdeObject object ---\n")
    validObject(object@data)
    validObject(object@sim.data)
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="NpdeObject",
  definition= function (.Object,data,sim.data,options=list(),prefs=list()){
    .Object@data<-data
    .Object@sim.data<-sim.data
    .Object@results<-new(Class="NpdeRes")
    .Object@results["ntot.obs"]<-data["ntot.obs"]
    .Object@results["not.miss"]<-data["not.miss"]
    .Object@results["icens"]<-data["icens"]
    opt<-npdeControl()
    if(length(options)>0) {
      for(i in names(options)) {
      	if(length(grep(i,names(opt)))==0) cat("Option",i, "not found, check spelling\n") else
      	opt[i]<-options[i]
      }
      i1<-grep("namsav",names(options))
      if(length(i1)!=0 && !is.na(i1)) {
        opt$namres<-paste(options[i1],".npde",sep="")
        opt$namgr<-paste(options[i1],".",opt$type.graph,sep="")
      }
    }
# Checking options
    opt<-check.control.options(opt)
    .Object@options<-opt
    graph.opt<-set.plotoptions(.Object)
    if(length(prefs)>0) {
    	for(i in names(prefs)) {
    		if(length(grep(i,names(graph.opt)))==0) cat("Graphical option",i, "not found, check spelling\n") else	graph.opt[i]<-prefs[i]
    	}
    }
    .Object@prefs<-graph.opt
    # Object validation
    validObject(.Object)
    return (.Object )
  }
)

###########################	Default options		#############################

#' Set options for an NpdeObject
#'
#' Set, replace and check options for an NpdeObject
#'
#' @name npdeControl
#' @aliases npdeControl replace.control.options check.control.options
#' @usage npdeControl(boolsave = TRUE, namsav = "output", type.graph = "eps", verbose = FALSE, calc.npde = TRUE, calc.pd = TRUE, decorr.method = "cholesky", cens.method = "omit", ties = TRUE, sample = FALSE)
#' @param boolsave whether to save the results (a file containing the numerical results and a file with the graphs)
#' @param namsav the root name of the files to save to (the file with the results will be named ROOTNAME.npde and the graphs will be saved to ROOTNAME.format where format is given by the type.graph argument)
#' @param type.graph type of graph to save to (one of "eps", "pdf", "jpeg", "png")
#' @param verbose a boolean; if TRUE, a message is printed as the computation of the npde begins for each new subject
#' @param calc.pd a boolean; TRUE to compute pd
#' @param calc.npde a boolean; TRUE to compute npde
#' @param decorr.method the method used to decorrelate simulated and observed data (see \code{\link{npde.decorr.method}})
#' @param cens.method the method used to handle censored data (see \code{\link{npde.cens.method}})
#' @param ties if FALSE, a smoothing will be applied to prediction discrepancies to avoid ties
#' @param sample if TRUE, the test on the pd will be performed after randomly sampling only pd per subject
#' @keywords methods

npdeControl<-function(boolsave=TRUE,namsav="output",type.graph="eps", verbose=FALSE,calc.npde=TRUE,calc.pd=TRUE,decorr.method="cholesky",cens.method="omit",ties=TRUE, sample=FALSE) {
# decorrelation methods:
#### cholesky: Cholesky decomposition
#### inverse: unique square root
#### polar: Cholesky decomposition combined with diagonalisation

# censoring methods: ECO TODO: find proper names
#### none: censored data removed, corresponding pd & npde set to NaN
#### pd.impute: when y<LOQ, sample pd in U(0,p_LOQ)
#### ipred: when y<LOQ, impute y as the model prediction and compute pd/npde for the completed dataset

# sample
#### when TRUE, for the tests based on pd, one sample per subject is randomly drawn 
   namres<-paste(namsav,".npde",sep="")
   namgr<-paste(namsav,".",type.graph,sep="")
   return(list(calc.pd=calc.pd,calc.npde=calc.npde,verbose=verbose, boolsave=boolsave,type.graph=type.graph,namsav=namsav,namres=namres,namgr=namgr, decorr.method=decorr.method,cens.method=cens.method,ties=ties,sample=sample))
}

replace.control.options<-function(opt,...) {
	args1<-match.call(expand.dots=TRUE)
	# These arguments are used by other functions and may be passed on via "..."
	legacy<-c("fix")
	if(length(args1)>2) {
		# General arguments: col, pch
		for(i in 3:length(args1)) {
			if(match(names(args1)[i],names(opt),nomatch=0)>0) {
				if(!is.null(eval(args1[[i]]))) opt[[names(args1)[i]]]<-eval(args1[[i]])
				} else {
					if(is.na(match(names(args1)[i],legacy))) cat("Argument",names(args1)[i],"not available, check spelling.\n")
				}
			}
		}
	opt<-check.control.options(opt)
	return(opt)
}

check.control.options<-function(opt) {
	if(!(opt$cens.method %in% c("omit","loq","ipred","ppred","fix","cdf"))) {
		cat("Warning: Censoring method",opt$cens.method,"unavailable, switching to default method (cdf)\n")
		opt$cens.method<-"cdf"
	}
	if(!(opt$decorr.method %in% c("cholesky","inverse","polar"))) {
		cat("Warning: Method",opt$decorr.method,"to decorrelate residuals is unavailable, switching to default method (cholesky)\n")
		opt$decorr.method<-"cholesky"
	}
	if(!(opt$type.graph %in% c("eps","png","pdf","jpeg"))) {
		cat("Warning: Type",opt$type.graph,"unrecognised; type of graph must be one of eps, png, pdf, jpeg, switching to default type (eps=Postcript)\n")
		opt$type.graph<-"eps"
	}
	for(bool.true in c("boolsave","verbose","calc.pd","calc.npde","ties")) {
		if(is.na(as.logical(opt[bool.true]))) {
			cat("Warning: Option",bool.true,"must be a logical (TRUE/FALSE), setting it to TRUE.\n")
			opt[bool.true]<-TRUE
		}
	}
	if(!is.logical(opt$sample)) {
		cat("Warning: Option",opt$sample,"must be a logical (TRUE/FALSE), setting it to FALSE.\n")
		opt$sample<-FALSE
	}
	invisible(opt)
}

##################################################################################

##' Get/set methods for NpdeData object
##' 
##' Access slots of a NpdeData using the object["slot"] format
##' 
##' @keywords methods
##' @exportMethod [

# Getteur
setMethod(
  f ="[",
  signature = "NpdeObject" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "data"={return(x@data)},
    "sim.data"={return(x@sim.data)},
    "results"={return(x@results)},
    "options"={return(x@options)},
    "prefs"={return(x@prefs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "NpdeObject" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "data"={x@data<-value},
    "sim.data"={x@sim.data<-value},
    "results"={x@results<-value},
    "options"={x@options<-value},
    "prefs"={x@prefs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
####				Summary method for NpdeObject			####
####################################################################################
setMethod("summary","NpdeObject",
  function(object) {
    cat("Object of class NpdeObject")
		if(length(object@data)>0) {
			cat(" containing the following main components\n")
			cat("   data: data\n")
    cat("       N=",object@data@N,"subjects\n")
    cat("       ntot.obs=",object@data@ntot.obs,"non-missing observations\n")
    cat("        subject id:",object@data@name.group,"\n")
    cat("        predictor (X):",object@data@name.predictor,"\n")
    cat("        response (Y):",object@data@name.response,"\n")
			if(length(object@data@name.covariates)>0)
				cat("        covariates:",object@data@name.covariates,"\n")
			if(length(object@sim.data)>0) {
				cat("   sim.data: simulated data: \n")
        cat("        number of replications: nrep=",object@sim.data@nrep,"\n")
			}
			if(length(object@results)>0) {
				cat("   results: results of the computation\n")
        cat("        ypred: individual predictions (E_i(f(theta_i,x)))\n")
        cat("        pd: prediction discrepancies\n")
				cat("        npde: normalised prediction distribution errors\n")
				cat("        ycomp: imputed responses for censored data\n")
				cat("        ploq: probability of being <LOQ for each observation\n")
			}
      cat("   options: options for computations\n")
      cat("   prefs: options for graphs\n")
		} else cat(", currently empty\n")
# The first elements of the list are the same as those returned previously, to maintain compatibility with the previous version of npde (1.2)
    obsdat<-data.frame(id=object@data@data[,object@data@name.group], xobs=object@data@data[,object@data@name.predictor], yobs=object@data@data[,object@data@name.response])
    addcol<-c(object@data@name.covariates,object@data@name.miss,object@data@name.cens,object@data@name.ipred)
    if(length(addcol)>0)
    	obsdat<-cbind(obsdat,object@data@data[,addcol])
    res<-list(obsdat=obsdat,ydobs=object@results@res$ydobs, ydsim=object@sim.data@datsim$ydsim,xerr=object@results@xerr, npde=object@results@res$npde,pd=object@results@res$pd, N=object@data@N, ntot.obs=object@data@ntot.obs, id=object@data@data[,object@data@name.group], x=object@data@data[,object@data@name.predictor], y=object@data@data[,object@data@name.response],nrep=object@sim.data@nrep, ypred=object@results@res$ypred,ycomp=object@results@res$ycomp, ploq=object@results@ploq,options=object@options,prefs=object@prefs)
    invisible(res)
 }
)

# Previous version of npde
# names(x)
# [1] "obsdat" "ydobs"  "ydsim"  "ypred"  "xerr"   "npde"   "pd"    

####################################################################################
####			Print and show methods for NpdeObject			####
####################################################################################
##################################################################################
# print/show/showall
# alias in class documentation

##' @S3method print NpdeObject

#setMethod("print","NpdeObject",
print.NpdeObject<-
function(x,nlines=10,...) {
    cat("Object of class NpdeObject\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    print(x@data,nlines=nlines)
    if(length(x@data@data)>0) {
      cat("\nSummary of original data:\n")
      cat("    vector of predictor",x@data@name.predictor,"\n")
      print(summary(x@data@data[,x@data@name.predictor]))
      cat("    vector of response",x@data@name.response,"\n")
      print(summary(x@data@data[,x@data@name.response]))
    }
    cat("-----------------------------------\n")
    cat("----         Key options       ----\n")
    cat("-----------------------------------\n")
    decmet<-c("Cholesky decomposition (upper triangular)","Inverse through diagonalisation","Polar method: Cholesky decomposition with diagonalisation")
    imet<-grep(x@options$decorr.method,c("cholesky","inverse","polar"))
    censmet<-c("Omit LOQ data (removed from data)", "Impute pd* and compute y* as F-1(pd*)","Impute y* to the LOQ", "Impute y* as the population model prediction", "Impute y* as the individual model prediction")
    icens<-grep(x@options$cens.method,c("omit","cdf","loq","ppred","ipred"))
    cat("Methods\n")
    cat("    compute prediction discrepancies (pd): ", ifelse(x@options$calc.pd,"yes","no"),"\n")
    cat("    compute normalised prediction distribution errors (npde): ", ifelse(x@options$calc.npde,"yes","no"),"\n")
    cat("    method for decorrelation: ",decmet[imet],"\n")
    cat("    method to treat censored data: ",censmet[icens],"\n")
    cat("Input/output\n")
    cat("    verbose (prints a message for each new subject): ", x@options$verbose,"\n")
    cat("    save the results to a file, save graphs:",x@options$boolsave,"\n")
    if(x@options$boolsave) {
      cat("    type of graph (eps=postscript, pdf=adobe PDF, jpeg, png):", x@options$type.graph,"\n")
      cat("    file where results should be saved: ",x@options$namres,"\n")
      cat("    file where graphs should be saved: ",x@options$namgr,"\n")
    }
    cat("-----------------------------------\n")
    if(length(x@results@res)>0) {
      cat("----      Results      ----\n")
      cat("-----------------------------------\n")
      print(x@results,nlines=nlines)
    } else cat("  No results\n")
  }
#)

setMethod("show","NpdeObject",
  function(object) {
#    cat("Object of class NpdeObject\n")
    cat("Object of class NpdeObject\n")
    cat("-----------------------------------------\n")
    cat("----        Component data           ----\n")
    cat("-----------------------------------------\n")
      show(object@data)
    cat("-----------------------------------------\n")
    if(length(object@results@res)>0) {    
      cat("----        Component results        ----\n")
      cat("-----------------------------------------\n")
      show(object@results)
    } else cat("  No results\n")
  }
)

##' @S3method showall NpdeObject
# Could be print, with only head of data
#setMethod("showall","NpdeObject",
showall.NpdeObject<-function(object) {
#    cat("Object of class NpdeObject\n")
    digits<-2;nsmall<-2
    cat("Object of class NpdeObject\n")
    cat("-----------------------------------\n")
    cat("----          Data             ----\n")
    cat("-----------------------------------\n")
    showall(object@data)
    cat("-----------------------------------\n")
    cat("----           Options         ----\n")
    cat("-----------------------------------\n")
    decmet<-c("Cholesky decomposition (upper triangular)","Inverse through diagonalisation")
    imet<-grep(object@options$decorr.method,c("cholesky","inverse"))
    censmet<-c("Omit LOQ data (removed from data)", "Impute pd* and compute y* as F-1(pd*)","Impute y* to the LOQ", "Impute y* as the population model prediction", "Impute y* as the individual model prediction")
    icens<-grep(object@options$cens.method,c("omit","cdf","loq","ppred","ipred"))
    cat("Methods\n")
    cat("    compute prediction discrepancies (pd): ", ifelse(object@options$calc.pd,"yes","no"),"\n")
    cat("    compute normalised prediction distribution errors (npde): ", ifelse(object@options$calc.npde,"yes","no"),"\n")
    cat("    method for decorrelation: ",decmet[imet],"\n")
    cat("    method to treat censored data: ",censmet[icens],"\n")
    cat("Input/output\n")
    cat("    verbose (prints a message for each new subject): ", object@options$verbose,"\n")
    cat("    save the results to a file, save graphs:",object@options$boolsave,"\n")
    if(object@options$boolsave) {
      cat("    type of graph (eps=postscript, pdf=adobe PDF, jpeg, png):", object@options$type.graph,"\n")
      cat("    file where results should be saved: ",object@options$namres,"\n")
      cat("    file where graphs should be saved: ",object@options$namgr,"\n")
    }
    cat("-----------------------------------\n")
    if(length(object@results@res)>0) {
      cat("----              Results            ----\n")
      cat("-----------------------------------------\n")
      showall(object@results)
    } else cat("  No results\n")
  }
#)

####################################################################################
####			subset method for NpdeObject				####
####################################################################################
##' @S3method subset NpdeObject

subset.NpdeObject<-function (x, subset, ...) {
    if (missing(subset)) 
        return(x)
    else {
        e <- substitute(subset)
	xdat<-as.data.frame(x["data"]["data"])
        r <- eval(e, xdat, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
# Data
    x1<-x
    x1["data"]["data"]<-x["data"]["data"][r,,drop=FALSE]
    if(length(x1["data"]["not.miss"])>0) {
      x1["data"]["not.miss"]<-x["data"]["not.miss"][r]
      x1["data"]["icens"]<-which(!x1["data"]["not.miss"])
    }
    id<-x1["data"]["data"][,x1["data"]["name.group"]]
    x1["data"]["N"]<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-c(nind.obs[match(unique(id),names(nind.obs))])
    x1["data"]["nind.obs"]<-nind.obs
    x1["data"]["ntot.obs"]<-length(id)
    x1["data"]["ind"]<-rep(1:x1["data"]["N"],times=nind.obs)
# Simulated data
    rsim<-rep(r,x["sim.data"]["nrep"])
    x1["sim.data"]["datsim"]<-x["sim.data"]["datsim"][rsim,,drop=FALSE] 
# Results
    x1["results"]["not.miss"]<-logical(0)
    x1["results"]["res"]<-x["results"]["res"][r,,drop=FALSE]
    x1["results"]["icens"]<-x1["data"]["icens"]
    x1["results"]["ntot.obs"]<-x1["data"]["ntot.obs"]
    if(length(x1["results"]["ploq"])>0) x1["results"]["ploq"]<-x["results"]["ploq"][r]
    if(length(x1["results"]["npde.sim"])>0) x1["results"]["npde.sim"]<-x["results"]["npde.sim"][r]
    if(length(x1["results"]["pd.sim"])>0) x1["results"]["pd.sim"]<-x["results"]["pd.sim"][r]
    return(x1)
}

####################################################################################
####				NpdeObject - test				####
####################################################################################

#' Test on npde or pd
#' 
#' Performs a global test on npde (default) or pd
#' 
#' @name gof.test
#' @aliases gof.test gof.test.numeric gof.test.NpdeRes gof.test.NpdeObject  print.gof.test
#' @param object an object (currently has methods for types numeric, NpdeRes and NpdeObject)
#' @param which whether the tests should be performed for npde (default), pd or npd (normalised pd) 
#' @param parametric whether parametric or non-parametric tests should be applied
#' @param \dots additional arguments passed on to the function; special arguments are \code{na.action}, which controls how to handle NAs in the results (\code{\link{na.action}}), \code{verbose} (if FALSE, suppresses printing of the results) and \code{covsplit} which requests the tests to be performed split by categories or quantiles of the data. If \code{covsplit} is TRUE, continuous covariates will be split in 3 categories (<Q1, Q1-Q3, >Q3) (see details in the PDF documentation), but this behaviour can be overriden by passing the argument \code{ncat=XXX} where XXX is the number of categories to divide the continuous covariates in.
#' @return A list with the following elements:
#' \describe{
#' \item{mean}{mean}
#' \item{se.mean}{standard error of the mean}
#' \item{var}{variance}
#' \item{se.var}{standard error on variance}
#' \item{kurtosis}{kurtosis (see \code{\link{kurtosis}})}
#' \item{skewness}{skewness (see \code{\link{skewness}})}
#' \item{p.value}{p-values for several tests (see below)}
#' }
#' @details If object is an NpdeObject and an argument covsplit=TRUE is given in \dots, in addition to the global descriptive statistics and tests, tests will be performed for each covariate in \code{which.cov}. This argument can be set in \dots; barring an explicit specification, the component \code{which.cov} of the prefs slot for a NpdeObject object will be used. The default value is \code{which.cov="all"}, which produces tests for each covariate in the dataset. Two additional dataframes will then be present:
#' \describe{
#' \item{cov.stat}{descriptive statistics and test p-values split by covariate and by categories}
#' \item{cov.p.value}{p-values split by covariate; for each covariate, two tests are performed: the first test is a correlation test for continuous covariates and a Chi-square test for categorical covariates; the second test is defined using the p-values of the global tests split by each category, and appling a Bonferroni correction to obtain an overall p-value (see PDF documentation for details)}
#' }
#' The p.value elements is a named vector with four components:
#' \describe{
#' \item{p.mean}{p-value for the mean test (Wilcoxon test if parametric=FALSE, Student test if parametric=TRUE)}
#' \item{p.var}{p-value for the variance test (parametric=FALSE, Fisher test if parametric=TRUE)}
#' \item{p.dist}{p-value for the distribution test (XXX if parametric=FALSE, XXX if parametric=TRUE)}
#' \item{p.global}{p-value for the global test (combination of the mean, variance and distribution tests with a Bonferroni correction)}
#' }
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
#' @references K. Brendel, E. Comets, C. Laffont, and F.Mentre. Evaluation of different tests based on observations for external model evaluation of  population analyses. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 37:49--65, 2010.
#' @seealso \code{\link{kurtosis}}, \code{\link{skewness}}
#' @keywords test
#' @examples
#' 
##' data(theopp)
##' 
##' @S3method gof.test NpdeObject
##' @export gof.test
##' @export gof.test.NpdeObject

gof.test.NpdeObject<-function(object,which="npde",parametric=TRUE, ...) {
	# Performs test on the selected variable (one of npde, pd or npd)
	if(length(object@results)==0) {
		cat("No results\n")
		return()
	}
	if(!which%in%c("pd","npde","npd")) {
		cat("Tests can be performed on one of: npde (default), pd, npd. Please choose one using the which argument.\n")
		return()
	}
	x<-switch(which,npde=object@results@res$npde,pd=object@results@res$pd, npd=qnorm(object@results@res$pd))
	x<-x[object["data"]["not.miss"]]
	args1<-match.call(expand.dots=TRUE)
	i1<-match("covsplit",names(args1))
	covsplit<-FALSE
	if(!is.na(i1) && !is.na(as.logical(as.character(args1[[i1]])))) covsplit<-as.logical(as.character(args1[[i1]]))
	i1<-match("sample",names(args1))
	sample.pd<-FALSE
	if(!is.na(i1) && !is.na(as.logical(as.character(args1[[i1]])))) sample.pd<-as.logical(as.character(args1[[i1]]))
	# 	verbose<-!(covsplit)
# 	i1<-match("verbose",names(args1))
# 	if(!is.na(i1) && !is.na(as.logical(args1[[i1]]))) verbose<-as.logical(args1[[i1]])
	if(covsplit & length(object["data"]["name.covariates"])>0) {
		which.cov<-object@prefs$which.cov
		i1<-match("which.cov",names(args1))
		if(!is.na(i1)) which.cov<-as.character(args1[[i1]])
		idx1<-!is.na(as.numeric(which.cov))
		if(length(idx1)>0) which.cov[idx1]<-object["data"]["name.covariates"][as.numeric(which.cov)[idx1]]
		if(which.cov[1]=="all") which.cov2<-object@data@name.covariates else
			which.cov2<-which.cov[match(which.cov, object@data@name.covariates,nomatch=0)>0]
		which.cov2<-unique(which.cov2)
		if(length(which.cov2)==0) {
			cat("No covariate matching names",which.cov,", performing a global test.\n")
			covsplit<-FALSE
		} else which.cov<-which.cov2
	} else covsplit<-FALSE
	if(sample.pd) {
		idobs<-object["data"]["data"][object["data"]["not.miss"], object["data"]["name.group"]]
		nind<-tapply(idobs,idobs,length) # individual numbers of observations (1xN)
		nind<-nind[match(unique(idobs),names(nind))]
		isamp<-cumsum(c(0,nind[-length(nind)]))+ceiling(runif(length(nind),0,nind))
		myres<-try(gof.test(x[isamp],which=which,parametric=parametric, ...))
		if(class(myres)!="try-error") print.gof.test(myres,which=as.character(which),...)
	} else {
		myres<-try(gof.test(x,which=which,parametric=parametric, ...))
		if(class(myres)!="try-error") print.gof.test(myres,which=as.character(which),...)
	}
	if(covsplit) {
		i1<-match("ncat",names(args1))
		ncat.cont<-NA
		if(!is.na(i1) && !is.na(as.integer(as.character(args1[[i1]])))) ncat.cont<-as.integer(as.character(args1[[i1]]))
		covtest<-glcov<-NULL
		for(icov in which.cov) {
			zecov<-object["data"]["data"][object["data"]["not.miss"],icov]
			idobs<-object["data"]["data"][object["data"]["not.miss"], object["data"]["name.group"]]
			ucov<-zecov[match(unique(idobs),idobs)]
#			if(!is.numeric(ucov) & length(unique(ucov))<=4) 
			if(is.factor(ucov)) {
# Categorical covariate: lm+anova if parametric, Kruskal-Wallis test otherwise
				covcont<-FALSE
				ncat<-length(unique(ucov))
				namcat<-paste(icov,sort(unique(zecov)),sep="=")
				if(parametric) {
					y<-try(anova(lm(x~zecov)))
					xval<-ifelse(class(y)=="try-error",NA,y$Pr[1])
					l1<-c(Covariate=icov,nb.categories=ncat,corr.pearson=xval)
				} else {
					y<-try(kruskal.test(x~zecov))
					xval<-ifelse(class(y)=="try-error",NA,y$p.value)
					l1<-c(Covariate=icov,nb.categories=ncat,corr.kruskal=xval)
				}
			} else {
# By default, test split by "<Q1","Q1-Q3",">Q3"
# If ncat.cont is given by the user, will be split by quantiles
				covcont<-TRUE
				if(is.na(ncat.cont)) ncat<-3 else ncat<-ncat.cont
				# Correlation test: Pearson if parametric, Spearman otherwise
				if(parametric) {
					y<-try(cor.test(zecov,x,method="pearson"))
					xval<-ifelse(class(y)=="try-error",NA,y$p.value)
					l1<-c(Covariate=icov,nb.categories=ncat,corr.pearson=xval)
				} else {
					y<-try(cor.test(zecov,x,method="spearman"))
					xval<-ifelse(class(y)=="try-error",NA,y$p.value)
					l1<-c(Covariate=icov,nb.categories=ncat,corr.spearman=xval)
				}
				if(is.na(ncat.cont)) {
					zecov<-cut(zecov,breaks=quantile(ucov,c(0,0.25,0.75,1)), include.lowest=TRUE, ordered_result=TRUE)
					namcat<-paste(icov,c("<Q1","Q1-Q3",">Q3"),sep=": ")
				} else {
						zecov<-cut(zecov,breaks=quantile(ucov,seq(0,1,length.out=(ncat+1))), include.lowest=TRUE, ordered_result=TRUE)
						namcat<-paste(icov,": Quantile ",0:(ncat-1),"/",ncat,"-",1:ncat,"/", ncat,sep="")
				}
			}
			icovtest<-NULL
			for(ic.cov in sort(unique(zecov))) {
				idx<-which(zecov==ic.cov & !is.na(x))
				if(sample.pd) {
					idx.obs<-idobs[idx]
					nind<-tapply(idx.obs,idx.obs,length) # individual numbers of observations (1xN)
					nind<-nind[match(unique(idx.obs),names(nind))]
					isamp<-cumsum(c(0,nind[-length(nind)]))+ceiling(runif(length(nind),0, nind))
					idx<-idx[isamp]
				}
				res1<-try(gof.test(x[idx], which=which,parametric=parametric, ...))
				if(class(res1)=="try-error") xval<-rep(NA,10) else xval<-unlist(res1)
				icovtest<-rbind(icovtest,c(ic.cov, namcat[which(sort(unique(zecov))==ic.cov)],xval))
			}
			for(i in 9:12) {
				xcal<-ncat*min(as.numeric(icovtest[,i]))
				l1<-c(l1,min(1,xcal))
			}
			glcov<-rbind(glcov,l1)
			covtest<-rbind(covtest,icovtest)
		}
		if(parametric) colnames(glcov)[4:7]<-c("  t-test                    ","  Fisher variance test      ","  SW test of normality      ", "Global adjusted p-value     ") else colnames(glcov)[4:7]<-c("  Wilcoxon signed rank test ","  Fisher variance test      ", "  SW test of normality      ","Global adjusted p-value     ")
		if(which=="pd") colnames(glcov)[6]<-"KS test of uniformity       "
		myres$cov.stat<-covtest
		myres$cov.p.value<-glcov
	}
	invisible(myres)
}

####################################################################################
####				NpdeObject - Main functions			####
####################################################################################

#' Main npde function
#' 
#' Main npde function, used to compute pd and npde
#' 
#' @param object A NpdeObject object
#' @return A NpdeObject object updated with the results
#' 
##' @keywords model internal
##' @S3method npde.main NpdeObject

#setMethod("npde.main","NpdeObject",
npde.main.NpdeObject<-
  function(object) {
#  	cat("Entering npde.main\n")
    if(object["options"]$calc.pd) object<-computepd(object)
    if(!object["options"]$calc.npde) gof.test(object,which="pd")
#    cat("Computation of pd successful\n")
    if(object["options"]$calc.npde) {
    	object<-computenpde(object)
    	gof.test(object)
    }
#    cat("Computation of npde successful\n")
    return(object)
  }
#)

#' Save the results contained in a NpdeObject object to a file
#' 
#' Save the results to a table on disk
#' 
##' @aliases npde.save,NpdeObject-method
#' @usage npde.save(object, ...)
#' @param object a NpdeObject object
#' @param \dots optional arguments to replace options in object
#' @details The following options can be changed by passing the appropriate arguments: namsav (string giving the root name of the files, an extension .npde will be added), nameres (string giving the full name of the file)
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
#' @keywords IO files
#' 
#' @export

setMethod("npde.save","NpdeObject",
  function(object, ...) {
  	args1<-match.call(expand.dots=TRUE)
  	i1<-match("namsav",names(args1))
  	i2<-match("namres",names(args1))
  	namres<-object["options"]$namres
  	if(!is.na(i2)) {
  		namres<-as.character(args1[[i2]])
  	} else {
  		if(!is.na(i1)) {
  		namres<-paste(as.character(args1[[i1]]),".npde",sep="")
  		}
  	}  	
  	if(namres=="") {
      cat("Please provide a filename in the namres item of the options component.\n")
      invisible()
    } else {
    	if(object@options$verbose) cat("Saving results in file",namres,"\n")
    }
    namcol<-c(object@data@name.group,object@data@name.predictor, object@data@name.response,object@data@name.cens,object@data@name.miss, object@data@name.covariates,object@data@name.ipred)
    saveres<-object["data"]["data"][,namcol]
    namcol2<-c("ypred","pd","npde")
    saveres<-cbind(saveres,object["results"]["res"][,namcol2])
    write.table(saveres,namres,row.names=FALSE,quote=FALSE)
  }
)

#' Save the graphs for a NpdeObject object to a file
#' 
#' Save the graphs to a file on disk
#' 
##' @aliases npde.graphs,NpdeObject-method
#' @usage npde.graphs(object, ...)
#' @param object a NpdeObject object
#' @param \dots optional arguments to replace options in object
#' @details The following options can be changed by passing the appropriate arguments: namsav (string giving the root name of the files, an extension depending on the type of graph will be added), namgr (string giving the full name of the file), type.graph (one of "eps", "pdf", "jpeg", "png")
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
#' @keywords IO files
#' 
#' @export

setMethod("npde.graphs","NpdeObject",
  function(object,...) {    
    args1<-match.call(expand.dots=TRUE)
    i1<-match("namgr",names(args1))
    i2<-match("namsav",names(args1))
    i3<-match("type.graph",names(args1))
    namsav<-object["options"]$namsav
    namgr<-object["options"]$namgr
    if(!is.na(i3)) {
      type.graph<-as.character(args1[[i3]])
    } else type.graph<-object["options"]$type.graph
    if(!is.na(i3) & is.na(i1) & is.na(i2)) namgr<-paste(namsav,type.graph,sep=".")
    if(!is.na(i1)) {
      if(is.na(i3)) namgr<-as.character(args1[[i1]]) else namgr<-as.character(args1[[i1]])
    } else {
      if(!is.na(i2)) {
        namsav<-as.character(args1[[i2]])
        namgr<-paste(namsav,type.graph,sep=".")
      }
    }
    if(length(namgr)==0 || namgr=="") {
      cat("Please provide a filename in the namgr item of the options component or as an option 'namgr=XXX' when calling npde.graphs().\n")
      return()
    } else {
    	if(object@options$verbose) cat("Saving graphs in file",namgr,"\n")
    }
    if(type.graph=="eps") postscript(namgr,onefile=TRUE,print.it=FALSE, horizontal=TRUE)
    if(type.graph=="jpeg") if(capabilities("jpeg")) jpeg(namgr) else {
         cat("R was not compiled with jpeg capabilities, switching to PDF format.\n")
	 type.graph<-"pdf"
	 namgr<-paste(namsav,type.graph,sep=".")
      }
    if(type.graph=="png") if(capabilities("png")) png(namgr) else {
         cat("R was not compiled with png capabilities, switching to PDF format.\n")
	 type.graph<-"pdf"
	 namgr<-paste(namsav,type.graph,sep=".")
      }
    if(type.graph=="pdf") pdf(namgr, onefile = TRUE)
    plot(object,...)
    dev.off()
  }
)


##################################################################################

#
#' Set graphical preferences
#'
#' This function is used to set options for graphs
#' 
#' @name set.plotoptions
##' @aliases set.plotoptions set.plotoptions,NpdeData-method set.plotoptions,NpdeObject-method set.plotoptions.NpdeData
#' @usage set.plotoptions(object, ...)
#' @param object an object of class NpdeData or NpdeObject
#' @param \dots arguments to replace default arguments (currently ignored)
#' @return a list of options for graphs
#' @details See documentation for a list of available options. 
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @keywords plot 

##' @S3method set.plotoptions NpdeObject
##' @export

set.plotoptions.NpdeObject<-function(object) {
	# setting default plot options
	plot.opt<-list(
		# General graphical options
		ask=FALSE,				# whether the program should ask before creating a new page
		new=TRUE,				# whether a new page should be created for the plot
		type.graph=object["options"]$type.graph,	# for compatibility with npde.graphs function
		interactive=FALSE,			# whether the user should be prompted before computing predictions or performing simulations for VPC, npde and wres
		main="",				# title
		sub="",					# sub-title
		xlab="",
		ylab="",
		cex=1,
		cex.axis=1,
		cex.lab=1,
		cex.main=1,
		mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
		xlim=c(),
		ylim=c(),
		xaxt=NULL, # A character which specifies the x axis type. Specifying "n" suppresses plotting.
		yaxt=NULL, # A character which specifies the y axis type. Specifying "n" suppresses plotting.
		axes=TRUE,				# Whether to plot the axes
		frame.plot=TRUE,	# Whether to add a box around the plotting region
		xlog=FALSE,
		ylog=FALSE,
		# Options for plot types
		ilist=c(1:object["data"]["N"]),
		plot.obs=TRUE,			# Whether observations, pd/ndpe should be plotted on top of the prediction bands (only applies if bands=TRUE)
		plot.loq=TRUE,			# Whether data under the LOQ should be plotted
		line.loq=TRUE,			# Whether an horizontal line should be plotted at Y=LOQ in data and VPC plots
		#		impute.loq=!(object["options"]$cens.method=="omit"),		# When TRUE, the imputed values are plotted for data under the LOQ; defaults to FALSE for "omit" method and TRUE for the other methods
		impute.loq=TRUE,		# When TRUE, the imputed values are plotted for data under the LOQ; defaults to TRUE
		smooth=FALSE,
		line.smooth="s",
		box=FALSE,				# Whether boxplots should be performed instead of scatterplots 
		which.cov="all",			# which covariates to plot 
		ncat=3,				# number of categories to bin continuous covariates
		which.resplot=c("res.vs.x","res.vs.pred","dist.qqplot","dist.hist"), # which type of residual plots
		# Colours, line types
		col="black",		# default colour of plot
		lty=1,
		lwd=1,
		type="b",
		pch=20,	 				# default symbol
		pch.pobs=20,				# symbol for observed/censored data
		pch.pcens=8,
		col.pobs="steelblue4",		# colour to plot observed/censored data
		col.pcens="red",
		col.lobs="steelblue4",
		lty.lobs=1,
		lwd.lobs=1,
		col.abline="DarkBlue",			# Characteristics of grid lines
		lty.abline=2,
		lwd.abline=2,
		boxwex=0.2,				# factor to scale width in boxplots
		varwidth=TRUE,			# use relative width for boxplots
		# Colours for prediction bands: VPC, npde, distribution plots
		range=3,
		bands=TRUE,				# Whether prediction bands should be added to the plots (including VPC)
		approx.pi=TRUE,			# Whether approximate prediction bands should be obtained for the distribution plots (see documentation)
		col.fillmed="pink",					# med: median
		col.fillpi="slategray1",		# pi: boundary of prediction interval
		col.lmed="indianred4",			# pop: description of observed data (eg median of observed concentrations)
		col.lpi="slategray4",
		lty.lmed=2,
		lty.lpi=2,
		lwd.lmed=1,
		lwd.lpi=1,
		# Binning options
		vpc.method="equal",			# method (one of "equal"=same nb of points in each interval, "width"=equally spaced intervals (on the log-scale if xlog=TRUE), "user"=user-defined breaks, "optimal"=Marc's optimal binning algorithm); for "user", the breaks must be specified in vpc.breaks (otherwise defaults back to "equal"), while for the other methods the number of bins must be specified in vpc.bin
		vpc.bin=10,				# nb of bins; the coordinates of the point used to summarise the data in each bin are the mean of the X and Y values of all points within the bins.
		vpc.breaks=NULL,			# user-defined breaks
		vpc.extreme=NULL, # can be set to a vector of 2 values to fine-tune the behaviour of the binning algorithm at the boundaries; specifying c(0.01,0.99) with the "equal" binning method and vpc.bin=10 will create 2 extreme bands containing 1% of the data on the X-interval, then divide the region within the two bands into the remaining 8 intervals each containing the same number of data; in this case the intervals will all be equal except for the two extreme intervals, the size of which is fixed by the user; complete fine-tuning can be obtained by setting the breaks with the vpc.method="user"
		vpc.interval=0.95, # size of the prediction intervals
		pi.size=0.95,			# width of the prediction interval on the quantiles
		vpc.beta=0.2,			# value of beta used to compute the variance-based criterion (Jopt,beta(I)) in the clustering algorithm
		vpc.lambda=0.3,			# value of lambda used in the penalised criterion to select the number of bins (if vpc.bin=NULL)
		bands.nrep=200)			# number of simulated datasets used to compute prediction bands
	
	plot.opt$name.X<-object["data"]["name.predictor"]
	plot.opt$name.Y<-object["data"]["name.response"]
	plot.opt$xlab<-paste(plot.opt$name.X," (",object["data"]["units"]$x,")", sep="")
	plot.opt$ylab<-paste(plot.opt$name.Y," (", object["data"]["units"]$y,")",sep="")
	return(plot.opt)
}

#
#' Plots a NpdeObject object
#'
#' Plots the data and diagnostic plots in a NpdeObject object
#' 
#' @param x a NpdeObject object
#' @param y unused, here for compatibility with the base plot function
#' @param \dots additional graphical parameters, which when given will supersede graphical preferences stored in the object
#' @details The default plot 
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.Mentre. Metrics for external model evaluation with an application to the population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49, 2006.
##' @seealso \code{\link{set.plotoptions}}
#' @keywords plot
#' @examples
#' 
#' data(theopp)
#' data(simtheopp)
#' 
#' x<-autonpde(theopp,simtheopp,iid="ID",ix="Time", iy="Conc", boolsave=FALSE)
#' plot(x)
#' 
##' @importFrom graphics plot
##' @method plot NpdeObject
#' @export

#	setMethod(f="plot",signature="NpdeObject", def=function(x,y,...) {
plot.NpdeObject <- function(x, y, ...) {
  	args1<-match.call(expand.dots=TRUE)
    i1<-match("new",names(args1))
    if(!is.na(i1)) force.new<-as.logical(as.character(args1[[i1]])) else force.new<-NULL
  	verbose<-x@options$verbose
  	i1<-match("verbose",names(args1))
  	if(!is.na(i1) && is.logical((as.character(args1[[i1]])))) verbose<-as.logical(as.character(args1[[i1]]))
  	i1<-match("plot.type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"default"
    if(verbose) cat("Selected plot type:",plot.type,"\n")
    i1<-match("which",names(args1))
    if(!is.na(i1)) {
      typmet<-as.character(args1[[i1]])
      typmet<-plot.type[plot.type!="c"]
    } else typmet<-"npde"
#    cat("plot.type=",plot.type,"\n")
    if(plot.type[1]=="default") {
#      cat("Creating default plots\n")
      plot.type<-c("qqplot","histogram","x.scatter", "pred.scatter")
      if(is.null(force.new) || force.new) par(mfrow=c(2,2),ask=x["prefs"]$ask)
      if(is.null(force.new)) force.new<-FALSE
    }
    pltyp<-c("data","ecdf","qqplot","histogram","x.scatter","pred.scatter", "cov.scatter","cov.x.scatter","cov.pred.scatter","cov.hist","cov.qqplot", "cov.ecdf","vpc","loq")
    ifnd<-pmatch(plot.type,pltyp)
    if(sum(is.na(ifnd))>0) {
      cat("The following plot types were not found or are ambiguous:", plot.type[is.na(ifnd)],"\n")
    }
    ifnd<-ifnd[!is.na(ifnd)]
    if(length(ifnd)==0) return("Plot type not found\n")
    plot.type<-pltyp[ifnd]
    interactive<-x["prefs"]$interactive
    namObj<-deparse(substitute(x))
# Check if pd or npde are present in the dataset, if not, perform the computation (ECO TODO remove ?)
    if(length(plot.type)>1 || plot.type[1]!="data") {
      icompute.pd<-icompute.npde<-icompute<-FALSE
      if(typmet %in% c("both","pd") & length(x["results"]["res"]$pd)==0) icompute.pd<-TRUE
      if(typmet %in% c("both","npde") & length(x["results"]["res"]$npde)==0) icompute.npde<-TRUE
      if(icompute.npde | icompute.pd) {
        icompute<-TRUE
        if(interactive) {
          i2<-(icompute.pd & icompute.npde)
          cok<-readline(prompt=paste("Computations will be performed to obtain", ifelse(icompute.pd,"pd",""),ifelse(i2," and ",""), ifelse(icompute.npde,"npde",""),", proceed ? (y/Y) [default=yes] ",sep=""))
          if(cok!="y"&cok!="Y"&cok!="yes"&cok!="") icompute<-FALSE
        } else 
        	{if(verbose) cat("Missing some elements for plots, will perform computations\n")}
      }
      if(icompute) {
      x["options"]$calc.pd<-icompute.pd
      x["options"]$calc.npde<-icompute.npde
      x<-npde.main(x)
      assign(namObj,x,envir=parent.frame())
      }
    }
#    cat(namObj,"\n")
  	if(typmet!="npde") x@prefs$bands<-FALSE
  	for(ipl in plot.type) {
      switch (EXPR=ipl,
    "data"={
    	if(verbose) cat("Plotting the data\n")
       npde.plot.data(x,...)
    },
    "x.scatter"={
    	if(verbose) cat("Plotting scatterplot versus X\n")
      npde.plot.scatter(x,xaxis="x",new=force.new,...)
    },
    "pred.scatter"={
    	if(verbose) cat("Plotting scatterplot versus predictions\n")
      npde.plot.scatter(x,xaxis="pred",new=force.new,...)
    },
    "cov.scatter"={
    	if(verbose) cat("Plotting scatterplot versus covariates\n")
    	npde.plot.scatter(x,xaxis="cov",new=force.new,...)
    },
    "qqplot"={
    	if(verbose) cat("Plotting QQ-plot of the distribution\n")
    	npde.plot.dist(x,dist.type="qqplot",new=force.new,...)
    },
    "histogram"={
    	if(verbose) cat("Plotting histogram of the distribution\n")
    	npde.plot.dist(x,dist.type="hist",new=force.new,...)
    },
    "ecdf"={
    	if(verbose) cat("Plotting the empirical distribution function of residuals\n")
    	npde.plot.dist(x,dist.type="ecdf",new=force.new,...)
    },
    "vpc"={
#      if(length(x["results"]["res"]$npde)>0) {
    	if(verbose) cat("Plotting VPC\n")
    	ich<-0
    	if(!x@prefs$bands) {x@prefs$bands<-TRUE;ich<-1}
      npde.plot.vpc(x,...)
    	if(ich==1) x@prefs$bands<-FALSE
#      }
    },
    "cov.x.scatter"={
    	if(verbose) cat("Plotting scatterplot versus X, splitted by covariate(s)\n")
      npde.plot.scatter(x,xaxis="x",covsplit=TRUE,...)
    },
    "cov.pred.scatter"={
    	if(verbose) cat("Plotting scatterplot versus predictions, splitted by covariate(s)\n")
      npde.plot.scatter(x,xaxis="pred",covsplit=TRUE,...)
    },
    "cov.hist"={
    	if(verbose) cat("Plotting histogram of the distribution, splitted by covariate(s)\n")
      npde.plot.dist(x,dist.type="hist",covsplit=TRUE,...)
    },
    "cov.qqplot"={
    	if(verbose) cat("Plotting histogram of the distribution, splitted by covariate(s)\n")
      npde.plot.dist(x,dist.type="qqplot",covsplit=TRUE,...)
    },
    "cov.ecdf"={
    	if(verbose) cat("Plotting histogram of the distribution, splitted by covariate(s)\n")
      npde.plot.dist(x,dist.type="ecdf",covsplit=TRUE,...)
    },
    "loq"={
      if(length(x["results"]["ploq"])>0) {
      	if(verbose) cat("Plotting p_LOQ=p(yobs<LOQ) \n")
        npde.plot.loq(x,...)
      }
    },
    cat("Plot ",ipl," not implemented yet\n")
     )
   }
}
#)

####################################################################################
