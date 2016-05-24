#' Capture-recapture model fitting function
#' 
#' Fits user specified models to some types of capture-recapture wholly in R
#' and not with MARK.  A single function that processes data, creates the
#' design data, makes the crm model and runs it
#' 
#' This function is operationally similar to the function \code{mark in RMark}
#' in that is is a shell that calls several other functions to perform the
#' following steps: 1) \code{\link{process.data}} to setup data and parameters
#' and package them into a list (processed data),2)
#' \code{\link{make.design.data}} to create the design data for each parameter
#' in the specified model, 3) \code{\link{create.dm}} to create the design
#' matrices for each parameter based on the formula provided for each
#' parameter, 4) call to the specific function for model fitting (now either
#' \code{\link{cjs}} or \code{\link{js}}). As with \code{mark} the calling
#' arguments for \code{crm} are a compilation of the calling arguments for each
#' of the functions it calls (with some arguments renamed to avoid conflicts).
#' If data is a processed dataframe (see \code{\link{process.data}}) then it
#' expects to find a value for \code{ddl}.  Likewise, if the data have not been
#' processed, then \code{ddl} should be NULL.  This dual calling structure
#' allows either a single call approach for each model or alternatively for the
#' data to be processed and the design data (\code{ddl}) to be created once and
#' then a whole series of models can be analyzed without repeating those steps.
#' 
#' There are some optional arguments that can be used to set initial values and
#' control other aspects of the optimization.  The optimization is done with
#' the R package/function \code{optimx} and the arguments \code{method} and
#' \code{hessian} are described with the help for that function.  In addition,
#' any arguments not matching those for \code{cjs} (the ...) are passed to
#' \code{optimx} allowing any of the other parameters to be set.  If you set
#' \code{debug=TRUE}, then at each function evaluation (\code{\link{cjs.lnl}}
#' the current values of the parameters and -2*log-likelihood value are output.
#' 
#' In the current implementation, a logit link is used to constrain the
#' parameters in the unit interval (0,1) except for probability of entry which
#' uses an mlogit and N which uses a log link. For the probitCJS model, a probit link is
#' used for the parameters. These could be generalized to
#' use other link functions. Following the notation of MARK, the parameters in
#' the link space are referred to as \code{beta} and those in the actual
#' parameter space of \code{Phi} and \code{p} as reals.
#' 
#' Initial values can be set in 2 ways.  To set a baseline intial value for the
#' intercept of \code{Phi} \code{p} set those arguments to some real value in
#' the open interval (0,1). All non-intercept beta parameters are set to zero.
#' Alternatively, you can specify in \code{initial}, a vector of initial values
#' for the beta parameters (on the logit scale).  This is most easily done by
#' passing the results from a previous model run using the result list element
#' \code{beta} as described below.  The code will match the names of the
#' current design matrix to the names in \code{beta} and use the appropriate
#' initial values. Any non-specified values are set to 0.  If there are no
#' names associated with the \code{initial} vector then they are simply used in
#' the specified order. If you do not specify initial values it is equivalent
#' to setting \code{Phi} and \code{p} to 0.5.
#' 
#' If you have a study with unequal time intervals between capture occasions,
#' then these can be specified with the argument \code{time.intervals}.
#' 
#' The argument \code{accumulate} defaults to \code{TRUE}.  When it is
#' \code{TRUE} it will accumulate common capture histories that also have
#' common design and common fixed values (see below) for the parameters.  This
#' will speed up the analysis because in the calculation of the likelihood
#' (\code{\link{cjs.lnl}} it loops over the unique values. In general the
#' default will be the best unless you have many capture histories and are
#' using many individual covariate(s) in the formula that would make each entry
#' unique. In that case there will be no effect of accumulation but the code
#' will still try to accumulate. In that particular case by setting
#' \code{accumulate=FALSE} you can skip the code run for accumulation.
#' 
#' Most of the arguments controlling the fitted model are contained in lists in
#' the arguments \code{model.parameters} and \code{design.parameters} which are
#' similar to their counterparts in \code{mark inb RMark}. Each is a named list
#' with the names being the parameters in the model (e.g., Phi and p in "cjs"
#' and "Phi","p","pent","N" in "js"). Each named element is also a list
#' containing various values defining the design data and model for the
#' parameter. The elements of \code{model.parameters} can include
#' \code{formula} which is an R formula to create the design matrix for the
#' parameter and \code{fixed} is a matrix of fixed values as described below.
#' The elements of \code{design.parameters} can include \code{time.varying},
#' \code{fields}, \code{time.bins},\code{age.bins}, and \code{cohort.bins}. See
#' \code{\link{create.dmdf}} for a description of the first 2 and
#' \code{\link{create.dm}} for a description of the last 3.
#' 
#' Real parameters can be set to fixed values using \code{fixed=x} where x is a
#' matrix with 3 columns and any number of rows.  The first column specifies
#' the particular animal (capture history) as the row number in the dataframe
#' \code{x}.  The second specifies the capture occasion number for the real
#' parameter to be fixed.  For \code{Phi} and \code{pent} these are 1 to
#' \code{nocc}-1 and for \code{p} they are 2 to \code{nocc} for "cjs" and 1 to
#' \code{nocc} for "js". This difference is due to the parameter labeling by
#' the beginning of the interval for Phi (e.g., survival from occasion 1 to 2)
#' and by the occasion for \code{p}.  For "cjs" \code{p} is not estimated for
#' occasion 1. The third element in the row is the real value in the closed
#' unit interval [0,1] for the fixed parameter.  This approach is completely
#' general allowing you to fix a particular real parameter for a specific
#' animal and occasion but it is a bit kludgy. Alternatively, you can set fixed
#' values by specifying values for a field called fix in the design data for a parameter.
#' If the value of fix is NA the parameter is estimated and if it is not NA then the real
#' parameter is fixed at that value.  If you also specify fixed as decribed above, they will over-ride any 
#' values you have also set with fix in the design data. To set all of the real values for a
#' particular occasion you can use the following example with the dipper data
#' as a template:
#' 
#' \code{model.parameters=list(Phi=list(formula=~1,}
#' \code{fixed=cbind(1:nrow(dipper),rep(2,nrow(dipper)),rep(1,nrow(dipper)))))}
#' 
#' The above sets \code{Phi} to 1 for the interval between occasions 2 and 3
#' for all animals. 
#' 
#' Alternatively, you could do as follows:
#' 
#' data(dipper)
#' dp=process.data(dipper)
#' ddl=make.design.data(dp)
#' ddl$Phi$fix=ifelse(ddl$Phi$time==2,1,NA)
#' 
#' At present there is no modification of the parameter count
#' to address fixing of real parameters except that if by fixing reals, a beta is not needed in the design it will be dropped.
#' For example, if you were to use ~time for Phi with survival fixed to 1 for time 2, then then beta for that time would not
#' be included.
#' 
#' To use ADMB (use.admb=TRUE), you need to install: 1) the R package R2admb, 2) ADMB, and 3) a C++ compiler (I recommend gcc compiler).
#' The following are instructions for installation with Windows. For other operating systems see (http://www.admb-project.org/downloads) and 
#'  (http://www.admb-project.org/tools/gcc/). 
#' 
#' Windows Instructions:
#'
#'  1) In R use install.packages function or choose Packages/Install Packages from menu and select R2admb.
#' 
#'  2) Install ADMB 11: http://admb-project.googlecode.com/files/admb-11-mingw-gcc4.5-32bit.exe. Put the software in C:/admb to
#'  avoid problems with spaces in directory name and for the function below to work.
#' 
#'  3) Install gcc 4.5 from: http://www.admb-project.org/tools/gcc/gcc452-win32.zip/view. Put in c:/MinGW
#' 
#' I use the following function in R to setup R2admb to access ADMB rather than adding to my path so gcc versions
#' with Rtools don't conflict. 
#' 
#' prepare_admb=function()
#' {
#'   Sys.setenv(PATH = paste("c:/admb/bin;c:admb/utilities;c:/MinGW/bin;", 
#'         Sys.getenv("PATH"), sep = ";"))
#'     Sys.setenv(ADMB_HOME = "c:/admb")
#'     invisible()
#' }
#' To use different locations you'll need to change the values used above
#' 
#' Before running crm with use.admb=T, execute the function prepare_admb().  You could put this function or the code it 
#' contains in your .First or .Rprofile so it runs each time you start R. 
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model (eg, "cjs", "js") 
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional vector of initial values for beta parameters; if
#' named from previous analysis only relevant values are used
#' @param groups Vector of names factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param method optimization method for function \code{optimx}
#' @param debug if TRUE, shows optimization output
#' @param hessian if TRUE, computes v-c matrix using hessian
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce
#' computation
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories and design matrices; amount used is 8*chunk_size/1e6 MB (default
#' 80MB)
#' @param control control string for optimization functions
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations for optimization 
#' @param scale vector of scale values for parameters
#' @param run if TRUE, it runs model; otherwise if FALSE can be used to test model build components 
#' @param burnin number of iterations for mcmc burnin; specified default not realistic for actual use
#' @param iter number of iterations after burnin for mcmc (not realistic default)
#' @param use.admb if TRUE creates data file for cjs.tpl and runs admb optimizer
#' @param crossed if TRUE it uses cjs.tpl or cjs_reml.tpl if reml=FALSE or TRUE respectively; if FALSE, then it uses cjsre which can use Gauss-Hermite integration
#' @param reml if TRUE uses restricted maximum likelihood
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param save.matrices for HMM models this option controls whether the gamma,dmat and delta matrices are saved in the model object
#' @param simplify if TRUE, design matrix is simplified to unique valus including fixed values
#' @param ... optional arguments passed to js or cjs and optimx
#' @importFrom graphics boxplot par
#' @importFrom stats as.formula binomial coef density
#'             glm.fit median model.frame model.matrix optim
#'              plogis pnorm predict rgamma rmultinom
#'              rnorm sd
#' @importFrom utils capture.output flush.console
#'             read.delim
#' @return crm model object with class=("crm",submodel) where submodel is
#' either "CJS" or "JS" at present.
#' @author Jeff Laake
#' @export crm
#' @import optimx Matrix Rcpp numDeriv
#' @useDynLib marked
#' @seealso \code{\link{cjs}}, \code{\link{js}},
#' \code{\link{make.design.data}},\code{\link{process.data}}
#' @keywords models
#' @examples
#' {
#' # cormack-jolly-seber model
#' # fit 3 cjs models with crm
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phit.pt=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)))
#' mod.Phit.pt
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,groups="sex",
#'    model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)))
#' mod.Phisex.pdot
#' ## if you have RMark installed you can use this code to run the same models 
#' ## by removing the comment symbol
#' #library(RMark)
#' #data(dipper)
#' #mod0=mark(dipper,
#' #model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),output=FALSE)
#' #summary(mod0,brief=TRUE)
#' #mod1=mark(dipper,
#' #model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)),output=FALSE)
#' #summary(mod1,brief=TRUE)
#' #mod2<-mark(dipper,groups="sex",
#' #model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)),output=FALSE)
#' #summary(mod2,brief=TRUE)
#' # jolly seber model
#' crm(dipper,model="js",groups="sex",
#'    model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)),accumulate=FALSE)
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # if you have RMark installed you can use this code to run the same models 
#' # by removing the comment 
#' #data(dipper)
#' #data(mstrata)
#' #mark(dipper,model.parameters=list(p=list(formula=~time)),output=FALSE)$results$beta
#' #mark(mstrata,model="Multistrata",model.parameters=list(p=list(formula=~1),
#' # S=list(formula=~1),Psi=list(formula=~-1+stratum:tostratum)),
#' # output=FALSE)$results$beta
#' #mod=mark(dipper,model="POPAN",groups="sex",
#' #   model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)))
#' #summary(mod)
#' #CJS example with hmm
#' crm(dipper,model="hmmCJS",model.parameters = list(p = list(formula = ~time)))
#' ##MSCJS example with hmm
#' data(mstrata)
#' ms=process.data(mstrata,model="hmmMSCJS",strata.labels=c("A","B","C"))
#' ms.ddl=make.design.data(ms)
#' ms.ddl$Psi$fix=NA
#' ms.ddl$Psi$fix[ms.ddl$Psi$stratum==ms.ddl$Psi$tostratum]=1
#' crm(ms,ms.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum)))
#' }
#' }
crm <- function(data,ddl=NULL,begin.time=1,model="CJS",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
 groups = NULL, time.intervals = NULL,debug=FALSE, method="BFGS", hessian=FALSE, accumulate=TRUE,chunk_size=1e7, 
 control=list(),refit=1,itnmax=5000,scale=NULL,run=TRUE,burnin=100,iter=1000,use.admb=FALSE,crossed=NULL,reml=FALSE,compile=FALSE,extra.args=NULL,
 strata.labels=NULL,clean=TRUE,save.matrices=TRUE,simplify=FALSE,...)
{
model=toupper(model)
ptm=proc.time()
if(is.null(crossed))crossed=FALSE
if(crossed)accumulate=FALSE
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
if(is.null(data$data))
{
   if(!is.null(ddl))
   {
      warning("Warning: specification of ddl ignored, as data have not been processed")
      ddl=NULL
   }
   message("Model: ",model,"\n")
   message("Processing data...\n")
   flush.console()
   data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=1, 
	   groups = groups, age.var = NULL, initial.ages = NULL, 
	   time.intervals = time.intervals,nocc=NULL,accumulate=accumulate,strata.labels=strata.labels)
}   
else
{
	data.proc=data
	model=data$model
}
#
# Setup parameter list
#
number.of.groups=1
if(!is.null(data.proc$group.covariates))number.of.groups=nrow(data.proc$group.covariates)
par.list=setup.parameters(data.proc$model,check=TRUE)
#
# Check validity of parameter list; stop if not valid
#
if(!valid.parameters(model,model.parameters)) stop()
parameters=setup.parameters(data.proc$model,model.parameters,data$nocc,number.of.groups=number.of.groups)
parameters=parameters[par.list]
# See if any formula contain random effects and set re
re=FALSE
for (i in 1:length(parameters))
{
	if(is.null(parameters[[i]]$formula)) parameters[[i]]$formula=~1
	mlist=proc.form(parameters[[i]]$formula)
	if(!is.null(mlist$re.model))re=TRUE
	if(parameters[[i]]$nointercept)parameters[[i]]$remove.intercept=TRUE
}
# currently if re, then set use.admb to TRUE
if(re) use.admb=TRUE
if(use.admb & !re) crossed=FALSE
# if re and accumulate=T, stop with message to use accumulate=FALSE
if(re & any(data.proc$freq>1)) stop("\n data cannot be accumulated (freq>1) with random effects; set accumulate=FALSE\n")
#
# If the design data have not been constructed, do so now
#
if(is.null(ddl)) 
{
	message("Creating design data...\n")
	flush.console()
	ddl=make.design.data(data.proc,design.parameters)
} else
{
	for (i in 1:length(parameters))
	{
		if(!is.null(ddl[[i]]$order))
		   if(any(ddl[[i]]$order!=1:nrow(ddl[[i]]))) 
			   stop(paste("Design data for parameter",names(parameters)[i],"is out of order."))
	}
    if(!is.null(design.parameters))
		for(parname in names(design.parameters))
			ddl$design.parameters[[parname]]=c(ddl$design.parameters[[parname]],design.parameters[[parname]])
	design.parameters=ddl$design.parameters
}
ddl=set.fixed(ddl,parameters) #   setup fixed values if old way used
if(simplify & !(substr(model,1,3)=="HMM"|(nchar(model)>=4 &substr(model,1,4)=="MVMS")))
{
	simplify=FALSE
	message("Can only use simplify with HMM models. simplify set to FALSE")
}
# Create design matrices for each parameter
dml=create.dml(ddl,model.parameters=parameters,design.parameters=design.parameters,chunk_size=chunk_size,simplify=simplify,use.admb=use.admb)
# For HMM call set.initial to get ptype and set initial values
if(substr(model,1,3)=="HMM"|(nchar(model)>=4 &substr(model,1,4)=="MVMS"))
	initial.list=set.initial(names(dml),dml,initial)
else
	initial.list=NULL
# if not running, return object with data,ddl,dml etc
if(!run) return(list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,ddl=ddl,dml=dml,results=initial.list))
# Depending on method set some values
if("SANN"%in%method)
{
	if(length(method)>1)
		warning("***SANN can only be used by itself; other methods ignored.")
	method="SANN"
    control$maxit=itnmax
}
if("nlminb"%in%method)
{
	control$eval.max=itnmax
	control$iter.max=itnmax
}
# Call estimation function which depends on the model
message("Fitting model\n")
if(model=="CJS")
    runmodel=cjs(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scale=scale,use.admb=use.admb,crossed=crossed,compile=compile,extra.args=extra.args,reml=reml,clean=clean,...)
if(model=="JS")
    runmodel=js(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=FALSE,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scale=scale,...)
if(model=="MSCJS")
	runmodel=mscjs(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				   refit=refit,control=control,itnmax=itnmax,scale=scale,use.admb=use.admb,re=re,compile=compile,extra.args=extra.args,clean=clean,...)
if(model=="PROBITCJS")
{
	if(is.null(initial))
	{
	    imat=process.ch(data.proc$data$ch,data.proc$data$freq,all=FALSE)
	    runmodel=probitCJS(ddl,dml,parameters=parameters,design.parameters=design.parameters,
		               imat=imat,iter=iter,burnin=burnin)
    }else
	    runmodel=probitCJS(ddl,dml,parameters=parameters,design.parameters=design.parameters,
					   initial=initial,iter=iter,burnin=burnin)	   
}
if(substr(model,1,3)=="HMM"|(nchar(model)>=4 &substr(model,1,4)=="MVMS"))
{
	if(substr(model,1,4)=="MVMS")
	{
		obslevels=data.proc$ObsLevels
		sup=data.proc$fct_sup(list(obslevels=obslevels))
	} else
		sup=NULL
	if(is.null(data.proc$strata.list) | substr(model,1,4)=="MVMS"){
		mx=data.proc$m
	}else{
		mx=list(ns=length(data.proc$strata.list$states),na=length(data.proc$strata.list[[names(data.proc$strata.list)[names(data.proc$strata.list)!="states"]]]))
	}
	runmodel=optimx(unlist(initial.list$par),HMMLikelihood,method=method,debug=debug,hessian=hessian,itnmax=itnmax,xx=data.proc$ehmat,mx=mx,
			        type=initial.list$ptype,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,control=control,
				    fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,
					parameters=parameters,sup=sup)
	par <- coef(runmodel, order="value")[1, ]
	runmodel=list(optim.details=as.list(summary(runmodel, order="value",par.select=FALSE)[1, ]))
	if(hessian)runmodel$hessian=attr(runmodel$optim.details,"details")$nhatend
	runmodel$convergence=runmodel$optim.details$convcode
	runmodel$options=list(accumulate=accumulate,initial=initial.list$par,method=method,
	                		chunk_size=chunk_size,itnmax=itnmax,control=control)
 	if(save.matrices)
	{
		runmodel$mat=HMMLikelihood(par=par,type=initial.list$ptype,xx=data.proc$ehmat,mx=mx,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
			fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,parameters=parameters,return.mat=TRUE,sup=sup)
	    if(model=="HMMCJS")
		{
			dimnames(runmodel$mat$gamma)[3:4]=list(c("Alive","Dead"),c("Alive","Dead"))
			dimnames(runmodel$mat$dmat)[3:4]=list(c("Missed","Seen"),c("Alive","Dead"))
		}else
		{
			dimnames(runmodel$mat$gamma)[3:4]=list(c(data.proc$strata.labels,"Dead"),c(data.proc$strata.labels,"Dead"))
			dimnames(runmodel$mat$dmat)[3:4]=list(data.proc$ObsLevels,c(data.proc$strata.labels,"Dead"))
		}
		names(dimnames(runmodel$mat$gamma))=c("Id","Occasion","From_state","To_state")
		names(dimnames(runmodel$mat$dmat))=c("Id","Occasion","Observation","State")
    }
	parlist=split(par,initial.list$ptype)
	par=vector("list",length=length(names(initial.list$par)))
	names(par)=names(initial.list$par)
	for(p in names(parlist))
	{
		par[[p]]=parlist[[p]]
		names(par[[p]])=colnames(dml[[p]]$fe)	
	}
	runmodel$beta=par
	runmodel$par=NULL
	runmodel$neg2lnl=2*runmodel$optim.details$value
	runmodel$AIC=runmodel$neg2lnl+2*sum(sapply(runmodel$beta,length))
	if(!is.null(runmodel$hessian))
	{
		runmodel$beta.vcv=solvecov(runmodel$hessian)$inv
		colnames(runmodel$beta.vcv)=names(unlist(runmodel$beta))
		rownames(runmodel$beta.vcv)=colnames(runmodel$beta.vcv)
	}
	class(runmodel)=c("crm","mle",model)
}
#
# Return fitted MARK model object or if external, return character string with same class and save file
if(!is.null(runmodel$convergence) && runmodel$convergence!=0&!use.admb)
{
	warning("******Model did not converge******")
	msg=attr(runmodel$optim.details,"details")$message
	if(is.null(msg)) msg="Exceeded maximum number of iterations"
	warning(msg)
}

object=list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,results=runmodel)
class(object)=class(runmodel)
if(!re & model!="MSCJS" & (nchar(model)<4 | (nchar(model)>=4 & substr(model,1,4)!="MVMS")))
   object$results$reals=predict(object,ddl=ddl,unique=TRUE,se=hessian)
cat(paste("\nElapsed time in minutes: ",round((proc.time()[3]-ptm[3])/60,digits=4),"\n"))
return(object)
}
# solvecov code was taken from package fpc: Christian
# Hennig chrish@@stats.ucl.ac.uk http://www.homepages.ucl.ac.uk/~ucakche/
solvecov=function (m, cmax = 1e+10)
# from package fpc
{
	options(show.error.messages = FALSE)
	covinv <- try(solve(m))
	if (class(covinv) != "try-error")
		coll = FALSE
	else {
		p <- nrow(m)
		cove <- eigen(m, symmetric = TRUE)
		coll <- TRUE
		if (min(cove$values) < 1/cmax) {
			covewi <- diag(p)
			for (i in 1:p) if (cove$values[i] < 1/cmax)
					covewi[i, i] <- cmax
				else covewi[i, i] <- 1/cove$values[i]
		}
		else covewi <- diag(1/cove$values, nrow = length(cove$values))
		covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
	}
	options(show.error.messages = TRUE)
	out <- list(inv = covinv, coll = coll)
	out
}

