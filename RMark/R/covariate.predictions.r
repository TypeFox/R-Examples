#' Compute estimates of real parameters for multiple covariate values
#' 
#' Computes real estimates for a dataframe of covariate values and the var-cov
#' matrix of the real estimates.
#' 
#' This function has a similar use as \code{\link{compute.real}} except that it
#' is specifically designed to compute real parameter estimates for multiple
#' covariate values for either a single model or to compute model averaged
#' estimates across a range of models within a marklist. This is particularly
#' useful for computing and plotting the real parameter as a function of the
#' covariate with pointwise confidence bands (see example below). The function
#' also computes a variance-covariance matrix for the real parameters.  For
#' example, assume you had a model with two age classes of young and adult and
#' survial for young was a function of weight and you wanted to estimate
#' survivorship to some adult age as a function of weight.  To do that you need
#' the survival for young as a function of weight, the adult survival, the
#' variance of each and their covariance.  This function will allow you to
#' accomplish tasks like these and many others.
#' 
#' When a variance-covariance matrix is computed for the real parameters, it
#' can get too large for available memory for a large set of real parameters.
#' Most models contain many possible real parameters to accomodate the general
#' structure even if there are very few unique ones. It is necessary to use the
#' most general structure to accomodate model averaging.  Most of the time you
#' will only want to compute the values of a limited set of real parameters but
#' possibly for a range of covariate values.  Use the argument \code{indices}
#' to select the real parameters to be computed.  The index is the value that
#' the real parameter has been assigned with the all-different PIM structure.
#' If you looked at the row numbers in the design data for the
#' \code{\link{dipper}} example, you would see that the parameter for p and Phi
#' are both numbered 1 to 21.  But to deal with multiple parameters effectively
#' they are given a unique number in a specific order.  For the CJS model, p
#' follows Phi, so for the dipper example, Phi are numbered 1 to 21 and then p
#' are numbered 22 to 42. You can use the function \code{\link{PIMS}} to lookup
#' the parameter numbers for a parameter type if you use
#' \code{simplified=FALSE}. For example, with the \code{\link{dipper}} data,
#' you could enter \code{PIMS(dipper.model,"p",simplified=FALSE)} and you would
#' see that they are numbered 22 to 42. Alternatively, you can use
#' \code{\link{summary.mark}} with the argument \code{se=TRUE} to see the list
#' of indices named \code{all.diff.index}.  They are included in a dataframe
#' for each model parameter which enables them to be selected based on the
#' attached data values (e.g., time, group etc).  For example, if you fitted a
#' model called \code{dipper.model} then you could use
#' \code{summary(dipper.model,se=TRUE)$real} to list the indices for each
#' parameter.
#' 
#' The argument \code{data} is a dataframe containing values for the covariates
#' used in the models.  The names for the fields should match the names of the
#' covariates used in the model.  If a time-varying covariate is used then you
#' need to specify the covariate name with the time index included as it is
#' specified in the data. You do not need to specify all covariates that were
#' used.  If a covariate in one or more models is not included in \code{data}
#' then the mean value will be used for each missing covariate.  That can be
#' very useful when you are only interested in prediction for one type of
#' parameters (eg Phi) when there are many covariates that are not interesting
#' in another parameter (e.g., p).  For each row in \code{data}, each parameter
#' specified in \code{indices} is computed with those covariate values.  So if
#' there were 5 rows in data and 10 parameters were specified there would be 10
#' sets of 5 (50) estimates produced.  If you do not want the full pairing of
#' data and estimates, create a field called \code{index} in \code{data} and
#' the estimate for that parameter will be computed with those specific values.
#' For example, if you wanted parameter 1 to be computed with 5 different
#' values of a covariate and then parameter 7 with 2 different covariate
#' values, you could create a dataframe with 5 rows each having an \code{index}
#' with the value 1 along with the relevant covariate data and an additional 2
#' rows with the \code{index} with the value 7.  If you include the field
#' \code{index} in \code{data} then you do not need to give a value for the
#' argument \code{indices}. However, if you are making the computations for parameters 
#' that use an mlogit link you must use the separate indices argument. If you try to
#' use the data.frame(index=...,cov) approach with mlogit parameters and you have
#' covariate values, the function will stop with an error.  Also, if you only include 
#' a portion of the indices in an mlogit set, it will also stop and issue an error and
#' tell you the set of indices that should be included for that mlogit set.  If you
#' were allowed to exclude some indices the result would be incorrect. 
#' 
#' @param model MARK model object or marklist
#' @param data dataframe with covariate values used for estimates; if it
#' contains a field called index the covariates in each row are only applied to
#' the parameter with that index and the argument indices is not needed; if data is not specified or
#' all individual covariate values are not specified, the mean individual covariate value is
#' used for prediction.
#' @param indices a vector of indices from the all-different PIM structure for
#' parameters to be computed (model.index value in the design data)
#' @param drop if TRUE, models with any non-positive variance for betas are
#' dropped
#' @param revised if TRUE it uses eq 6.12 from Burnham and Anderson (2002) for
#' model averaged se; otherwise it uses eq 4.9
#' @param mata if TRUE, create model averaged tail area confidence intervals as described by Turek and Fletcher
#' @param alpha The desired lower and upper error rate.  Specifying alpha=0.025
#' corresponds to a 95% MATA-Wald confidence interval, an' 
#' alpha=0.05 to a 90% interval.  'alpha' must be between 0 and 0.5.
#' Default value is alpha=0.025.
#' @param normal.lm Specify normal.lm=TRUE for the normal linear model case, and 
#' normal.lm=FALSE otherwise.  When normal.lm=TRUE, the argument 
#' 'residual.dfs' must also be supplied.  See USAGE section, 
#' and Turek and Fletcher (2012) for additional details.
#' @param residual.dfs A vector containing the residual (error) degrees of freedom 
#' under each candidate model.  This argument must be provided 
#' when the argument normal.lm=TRUE.
#' @param ... additional arguments passed to specific functions
#' @return A list is returned containing a dataframe of estimates, a
#' var-cov matrix, and a reals list: 
#' \item{estimates}{ data frame containing estimates, se,
#' confidence interval and the data values used to compute the estimates}
#' \item{vcv}{variance-covariance matrix of real estimates}
#' \item{reaks}{list of dataframes with the estimate and se used for each model}
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{compute.real}},\code{\link{model.average}}
#' @keywords utility
#' @examples
#' pathtodata=paste(path.package("RMark"),"extdata",sep="/")
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' #
#' # indcov1.R
#' #
#' # CJS analysis of the individual covariate data from 12.2 of
#' # Cooch and White
#' #
#' # Import data (indcov1.inp) and convert it from the MARK inp file 
#' # format to the RMark format using the function convert.inp  
#' # It is defined with 1 group but 2 individual covariates of mass and 
#' # squared.mass
#' #
#' indcov1=convert.inp(paste(pathtodata,"indcov1",sep="/"),
#'                covariates=c("mass","squared.mass"))
#' #
#' # Next create the processed dataframe and the design data.
#' #
#'   ind1.process=process.data(indcov1,model="CJS")
#'   ind1.ddl=make.design.data(ind1.process)
#' #
#' # Next create the function that defines and runs the set of models
#' # and returns a marklist with the results and a model.table.  
#' # It does not have any arguments but does use the ind1.process 
#' # and ind1.ddl objects created above in the workspace. The function 
#' # create.model.list is the function that creates a dataframe of the
#' # names of the parameter specifications for each parameter in that 
#' # type of model. If none are given for any parameter, the default
#' # specification will be used for that parameter in mark.  The 
#' # first argument of mark.wrapper is the model list of parameter 
#' # specifications.  Remaining arguments that are passed to
#' # mark must be specified using the argument=value specification 
#' # because the arguments of mark were not repeated in mark.wrapper 
#' # so they must be passed using the argument=value syntax.
#' #
#' ind1.models=function()
#' {
#'   Phi.dot=list(formula=~1)
#'   Phi.mass=list(formula=~mass)
#'   Phi.mass.plus.mass.squared=list(formula=~mass + squared.mass)
#'   p.dot=list(formula=~1)
#'   cml=create.model.list("CJS")
#'   results=mark.wrapper(cml,data=ind1.process,ddl=ind1.ddl,adjust=FALSE)
#'   return(results)
#' }
#' #
#' # Next run the function to create the models and store the results in
#' # ind1.results which is a marklist. Note that beta estimates will differ
#' # from Cooch and White results because we are using covariate values 
#' # directly rather than standardized values.
#' #
#' ind1.results=ind1.models()
#' #
#' # Next compute real parameter values for survival as a function of 
#' # mass which are model-averaged over the fitted models.
#' #
#' minmass=min(indcov1$mass)
#' maxmass=max(indcov1$mass)
#' mass.values=minmass+(0:30)*(maxmass-minmass)/30
#' Phibymass=covariate.predictions(ind1.results,
#'    data=data.frame(mass=mass.values,squared.mass=mass.values^2),
#'     indices=c(1))
#' #
#' # Plot predicted model averaged estimates by weight with pointwise 
#' # confidence intervals
#' #
#' plot(Phibymass$estimates$mass, Phibymass$estimates$estimate,
#'    type="l",lwd=2,xlab="Mass(kg)",ylab="Survival",ylim=c(0,.65))
#' lines(Phibymass$estimates$mass, Phibymass$estimates$lcl,lty=2)
#' lines(Phibymass$estimates$mass, Phibymass$estimates$ucl,lty=2)
#'
#' # indcov2.R
#' #
#' # CJS analysis of the individual covariate data from 12.3 of 
#' # Cooch and White
#' #
#' # Import data (indcov2.inp) and convert it from the MARK inp file 
#' # format to the RMark format using the function convert.inp  It is
#' # defined with 1 group but 2 individual covariates of mass and 
#' # squared.mass
#' #
#' indcov2=convert.inp(paste(pathtodata,"indcov2",sep="/"),
#'         covariates=c("mass","squared.mass"))
#' #
#' # Standardize covariates
#' #
#' actual.mass=indcov2$mass
#' standardize=function(x,z=NULL)
#' {
#'    if(is.null(z))
#'    {
#'        return((x-mean(x))/sqrt(var(x)))
#'    }else 
#'    {
#'       return((x-mean(z))/sqrt(var(z)))
#'    }
#' }
#' indcov2$mass=standardize(indcov2$mass)
#' indcov2$squared.mass=standardize(indcov2$squared.mass)
#' #
#' # Next create the processed dataframe and the design data.
#' #
#'   ind2.process=process.data(indcov2,model="CJS")
#'   ind2.ddl=make.design.data(ind2.process)
#' #
#' # Next create the function that defines and runs the set of models and 
#' # returns a marklist with the results and a model.table.  It does not
#' # have any arguments but does use the ind1.process and ind1.ddl 
#' # objects created above in the workspace. The function create.model.list 
#' # is the function that creates a dataframe of the names of the parameter
#' # specifications for each parameter in that type of model. If none are 
#' # given for any parameter, the default specification will be used for
#' # that parameter in mark.  The first argument of mark.wrapper is the
#' # model list of parameter specifications.  Remaining arguments that are
#' # passed to mark must be specified using the argument=value specification 
#' # because the arguments of mark were not repeated in mark.wrapper so 
#' # they must be passed using the argument=value syntax.
#' #
#' ind2.models=function()
#' {
#'   Phi.dot=list(formula=~1)
#'   Phi.time=list(formula=~time)
#'   Phi.mass=list(formula=~mass)
#'   Phi.mass.plus.mass.squared=list(formula=~mass + squared.mass)
#'   Phi.time.x.mass.plus.mass.squared=
#'         list(formula=~time:mass + time:squared.mass)
#'   Phi.time.mass.plus.mass.squared=
#'      list(formula=~time*mass + squared.mass+ time:squared.mass)
#'   p.dot=list(formula=~1)
#'   cml=create.model.list("CJS")
#'   results=mark.wrapper(cml,data=ind2.process,ddl=ind2.ddl,adjust=FALSE,threads=2)
#'   return(results)
#' }
#' #
#' # Next run the function to create the models and store the results in
#' # ind2.results which is a marklist. Note that beta estimates will differ
#' # because we are using covariate values directly rather than 
#' # standardized values.
#' #
#' ind2.results=ind2.models()
#' #
#' # Next compute real parameter values for survival as a function of 
#' # mass which are model-averaged over the fitted models. They are 
#' # standardized individually so the values have to be chosen differently.
#' #
#' minmass=min(actual.mass)
#' maxmass=max(actual.mass)
#' mass.values=minmass+(0:30)*(maxmass-minmass)/30
#' squared.mass.values=mass.values^2
#' mass.values=standardize(mass.values,actual.mass)
#' squared.mass.values=standardize(squared.mass.values,actual.mass^2)
#' Phibymass=covariate.predictions(ind2.results,
#'  data=data.frame(mass=mass.values,squared.mass=squared.mass.values),
#'   indices=c(1:7))
#' #
#' # Plot predicted model averaged estimates by weight with pointwise 
#' # confidence intervals
#' #
#' par(mfrow=c(4,2))
#' for (i in 1:7)
#' {
#'    mass=minmass+(0:30)*(maxmass-minmass)/30
#'    x=Phibymass$estimates
#'    plot(mass,x$estimate[x$par.index==i],type="l",lwd=2,
#'     xlab="Mass(kg)",ylab="Survival",ylim=c(0,1),main=paste("Time",i))
#'    lines(mass, x$lcl[x$par.index==i],lty=2)
#'    lines(mass, x$ucl[x$par.index==i],lty=2)
#' }
#' }
#' 
covariate.predictions <- function(model,data=NULL,indices=NULL,drop=TRUE, revised=TRUE, mata=FALSE, normal.lm=FALSE, residual.dfs=0, alpha=0.025,...)
{
# ------------------------------------------------------------------------------------------------
#
# covariate.predictions
#            -   computes real estimates and var-cov for a set of
#                covariate values
#
# Arguments:
#
#   model          - MARK model object or marklist
#   data           - dataframe with covariate values used for estimates
#                     if it contains a field index the covariates in each row
#                       are only applied to that index and the argument indices is not needed
#   indices        - all-different PIM indices for parameters to be calculated
#   drop           - if TRUE, models with any non-positive variance for betas are dropped
#
#
# Value (list):
#
#   estimates        - data frame containing estimates and se
#   vcv              - variance-covariance matrix of real estimates
#
# ------------------------------------------------------------------------------------------------
   if(class(model)[1]!="marklist")
   {
     number.of.models=1
     model=load.model(model)
     model.list=list(model)
   }
   else
   {
     number.of.models=length(model)-1
     if(is.null(model$model.table))
        stop("\nmarklist created by collect.models must contain a model.table to use model.average\n")
     model.list=model
     model.table=model$model.table
   }
#   if(is.null(data))
#      stop("\n data argument must be specified\n")
   if(!is.null(data)&!is.data.frame(data))
      stop("\n data argument must be a dataframe. Do not use processed data list.\n")
#
#   If there is an index field in data, then only use that row of data for that index
#   That filtering is done below after the complete design matrix is created for each model.
#
   if(!is.null(data$index))
   {
      index=data$index
      indices=index
      replicate.values=FALSE
   }
   else {
	   if(!is.null(data))
		   replicate.values=TRUE
	   else{
		   replicate.values=FALSE
		   if(is.null(indices))
			   stop("\nValue for indices argument must be given or index field must be included in data argument\n")
		   else
		   {
			   data=data.frame(index=indices)
			   index=indices
		   }
	   }
   }
#
# Determine if any of the models should be dropped because beta.var non-positive
#
if(number.of.models>1)
{
dropped.models=NULL
if(drop)
{
   for (i in 1:dim(model.table)[1])
   {
      model=load.model(model.list[[i]])
      model.indices=unique(model$simplify$pim.translation[indices])
      used.beta=which(apply(model$design.matrix[model.indices,,drop=FALSE],2,function(x)!all(x=="0")))
      if(any(diag(model$results$beta.vcv[used.beta,used.beta])<0))
      {
         dropped.models=c(dropped.models,i)
         message("\nModel ",i,"dropped from the model averaging because one or more beta variances are not positive\n")
      }
   }
#
# If any models have been dropped, recompute the weights for the model table
#
   if(!is.null(dropped.models))
   {
      model.table$weight[as.numeric(row.names(model.table))%in%dropped.models]=0
      model.table$weight=model.table$weight/sum(model.table$weight)
   }
}
}
else
{
	model.indices=unique(model$simplify$pim.translation[indices])
	used.beta=which(apply(model$design.matrix[model.indices,,drop=FALSE],2,function(x)!all(x=="0")))
	if(any(diag(model$results$beta.vcv[used.beta,used.beta])<0))
		message("\nModel has one or more beta variances that are not positive\n")
}
dropped.models=NULL
reals=vector("list",length=number.of.models)
firstmodel=TRUE
for (j in 1:number.of.models)
{
   if(j %in% dropped.models) next
#
#   Assign value of chat - overdispersion constant
#
   model=load.model(model.list[[j]])
   if(is.null(model$chat))
     chat=1
   else
     chat=model$chat
   beta=model$results$beta$estimate
   model.indices=unique(model$simplify$pim.translation[indices])
   links=NULL
   design=NULL
   fixedparms=NULL
   boundaryparms=NULL
#   
#  If there are no individual covariates in the model, the code below will extract only those
#  rows in the DM.  
#
   nmlogit=0
   if(is.null(model$covariates)&ncol(data)==1 && names(data)=="index")
#   if(ncol(data)==1 && names(data)=="index")
   {
     for (i in 1:nrow(data))
     {
        model.index=model$simplify$pim.translation[index[i]]   
        design=rbind(design,model$design.matrix[model.index,,drop=FALSE])
        if(length(model$links)==1)
           links=c(links,model$links)
        else
		   links=c(links,model$links[model.index])	
        if(!is.null(model$fixed))
        {
           if(is.null(model$simplify))
              xfixedparms=(1:dim(model$design.matrix)[1])%in%model$fixed$index
           else                                                           
              xfixedparms=(1:dim(model$design.matrix)[1])%in%model$simplify$pim.translation[model$fixed$index]
        }
        else
           xfixedparms=rep(FALSE,dim(model$design.matrix)[1])
        fixedparms=c(fixedparms,xfixedparms[model.index])
        boundaryparms=c(boundaryparms,(model$results$real$se==0 & !xfixedparms)[model.index])
     }
	 mlinks=substr(links,1,6)=="mlogit"
	 if(any(mlinks))
	 {
		 mlogit.list=split(1:length(links[mlinks]),links[mlinks])
		 for(i in 1:length(mlogit.list))
		 {
			 if(any(!which(model$links==names(mlogit.list)[i])%in%model.indices)) stop(paste("\n some indices are missing for parameters with",
								 names(mlogit.list)[i],"\n should include",paste(which(model$simplify$pim.translation%in%which(model$links==names(mlogit.list)[i])),collapse=",")))
		 }
	 }
   }
   else
# If there are data values for covariates then it fills in a complete DM for each record in the data
# and appends to create a large DM with all individual DMs stacked on each other.  Then in the 
# next loop if an index is specified in data, it selects out only the rows of the large DM for those 
# specific parameter indices.
#
   {
     for (i in 1:nrow(data))
     {
        xdesign=fill.covariates(model,find.covariates(model,data[i,,drop=FALSE]))[model.indices,,drop=FALSE]
        if(length(model$links)==1)
           links=c(links,rep(model$links,dim(xdesign)[1]))
        else
		if(length(model$links)==1)
			links=c(links,model$links)
		else
		{
			clinks=model$links[model.indices]
			mlinks=substr(clinks,1,6)=="mlogit"
			if(any(mlinks))
			{
				if(!replicate.values) 
					stop("Computations for mlogit parameters with covariate values cannot be specified with \nindex column in data; use separate indices argument" )
				mlogit.list=split(1:length(clinks[mlinks]),clinks[mlinks])
				for(i in 1:length(mlogit.list))
				{
					if(any(!which(model$links==names(mlogit.list)[i])%in%model.indices)) stop(paste("\n some indices are missing for parameters with",
										names(mlogit.list)[i],"\n should include",paste(which(model$simplify$pim.translation%in%which(model$links==names(mlogit.list)[i])),collapse=",")))
					nmlogit=nmlogit+1
					clinks[mlinks][mlogit.list[[i]]]=paste("mlogit(",nmlogit,")",sep="")
				}
			}	
			links=c(links,clinks)		
		}
		design=rbind(design,xdesign)
        if(!is.null(model$fixed))
        {
           if(is.null(model$simplify))
              xfixedparms=(1:dim(model$design.matrix)[1])%in%model$fixed$index
           else                                                           
              xfixedparms=(1:dim(model$design.matrix)[1])%in%model$simplify$pim.translation[model$fixed$index]
        }
        else
           xfixedparms=rep(FALSE,dim(model$design.matrix)[1])
        fixedparms=c(fixedparms,xfixedparms[model.indices])
        boundaryparms=c(boundaryparms,(model$results$real$se==0 & !xfixedparms)[model.indices])
     }
#
#    Filter rows if each covariate set is not replicated for each parameter
#
     if(!replicate.values)
     {
        row.numbers=NULL
        for(i in 1:length(index))
        {
          model.index=model$simplify$pim.translation[index[i]]
          row.numbers=c(row.numbers,match(model.index,model.indices)+(i-1)*length(model.indices))
        }
        design=design[row.numbers,,drop=FALSE]
        fixedparms=fixedparms[row.numbers]
        boundaryparms=boundaryparms[row.numbers]
        links=links[row.numbers]
     }
   }
#                                                                         
#   The following shouldn't happen unless the model and design matrices are mixed between models
#
   if(dim(design)[2]!=length(beta))
       stop("Mismatch between number of design columns and length of beta")
#
#  Convert to a numeric matrix
#
   design=matrix(as.numeric(design),nrow=dim(design)[1])
# 
#  Create vector of fixed parameters to cope with partial fixing of mlogits
#
   fixedvalues=rep(NA,nrow(design))
   fixedvalues[fixedparms]=model$fixed$value[match(indices[fixedparms],model$fixed$index)]
#
#  Compute real parameters; if neither se or vcv then return vector of real parameters
#
   real=convert.link.to.real(design%*%beta,links=links,fixed=fixedvalues)
#
#  Set fixed real parameters to their fixed values
#
   real[fixedparms]=fixedvalues[fixedparms]
#   Previously read as follows -- fixed in v1.9.5
#   real[fixedparms]=model$fixed$value[model$fixed$index%in%indices[fixedparms]]
#
#  Compute vc matrix for real parameters which needs to be modified for
#  any mlogit parameters
#
   if(length(model$links)==1)
      deriv.real=deriv_inverse.link(real,design,model$links)
   else
      deriv.real=t(apply(data.frame(real=real,x=design,links=links),1,
            function(x){deriv_inverse.link(as.numeric(x[1]),as.numeric(x[2:(length(x)-1)]),x[length(x)])}))
   vcv.real=deriv.real%*%model$results$beta.vcv%*%t(deriv.real)
#
#   If vcv=TRUE, compute v-c matrix and std errors of real estimates
#   To handle any mlogit parameters compute pseudo-real estimates using log in place of mlogit
#
   ind=grep("mlogit",links,ignore.case=TRUE)
   templinks=links
   if(length(ind)>0)
   {
      templinks[ind]="log"
      pseudo.real=as.vector(convert.link.to.real(design%*%beta,links=templinks))
      pseudo.real[fixedparms]=fixedvalues[fixedparms]
      pseudo.real[ind][fixedparms[ind]]=exp(pseudo.real[ind][fixedparms[ind]])
#
#     Compute first derivatives of pseudo-real (for any mlogit parameter)
#     estimates with respect to beta parameters
#
      if(length(templinks)==1)
        deriv.pseudo=deriv_inverse.link(pseudo.real,design,templinks)
      else
         deriv.pseudo=t(apply(data.frame(real=pseudo.real,x=design,links=templinks),1,
               function(x){deriv_inverse.link(as.numeric(x[1]),as.numeric(x[2:(length(x)-1)]),x[length(x)])}))
      deriv.pseudo[fixedparms,]=0
      vcv.pseudo=chat*deriv.pseudo%*%model$results$beta.vcv%*%t(deriv.pseudo)
#
#     Apply chain rule to get variance of real parameters which has mlogits
#     expressed as zi/(1+z1+...zk) where k is number of mlogit components-1 and
#     non-mlogits are expressed as zi.
#     bottom is either 1 for non-mlogits and the sum for mlogits
#     pbottom is partial with respect to zi
#
      mlogits=outer(links,links,function(x,y)as.numeric(x==y))*as.numeric(substr(links,1,6)=="mlogit" |substr(links,1,6)=="MLogit" )
      pbottom=matrix(0,nrow=dim(vcv.pseudo)[1],ncol=dim(vcv.pseudo)[1]) + mlogits
      bottom=diag(nrow=dim(vcv.pseudo)[1])*(1-as.numeric(substr(links,1,6)=="mlogit" |substr(links,1,6)=="MLogit"))+
          mlogits + pbottom*apply(pbottom*pseudo.real,2,sum)
      deriv.pseudo=(diag(nrow=dim(vcv.pseudo)[1])*bottom-pseudo.real*pbottom)/bottom^2
      deriv.pseudo[is.nan(deriv.pseudo)]=0
      vcv.real=deriv.pseudo%*%vcv.pseudo%*%t(deriv.pseudo)
   }
   else
      vcv.real=chat*vcv.real
#
#  Compute conf interval taking into account use of logit transform for mlogit links
#
   link.se=suppressWarnings(sqrt(chat*diag(design%*%model$results$beta.vcv%*%t(design))))
   link.se[is.na(link.se)]=0
   ind=unique(c(grep("mlogit",links,ignore.case=TRUE),which(links%in%c("sin","Sin","LogLog","loglog","CLogLog","cloglog"))))
   linkse=suppressWarnings(sqrt(diag(vcv.real)[ind])/(real[ind]*(1-real[ind])))
   linkse[is.na(linkse)]=0
   linkse[is.infinite(linkse)]=0
   link.se[ind]=linkse
   link.values=design%*%beta
   link.values[ind]=suppressWarnings(log(real[ind]/(1-real[ind])))
   link.values[ind][abs(real[ind]-1)<1e-7]=100
   link.values[ind][abs(real[ind]-0)<1e-7]=-100
   links[ind]="logit"
   real.lcl=convert.link.to.real(link.values-1.96*link.se,links=links)
   real.ucl=convert.link.to.real(link.values+1.96*link.se,links=links)
   real.lcl[fixedparms]=real[fixedparms]
   real.ucl[fixedparms]=real[fixedparms]   
#
#  Set v-c values of fixed parameters to 0
#
   vcv.real[fixedparms,]=0
   vcv.real[,fixedparms]=0
   se.real=suppressWarnings(sqrt(diag(vcv.real)))
   se.real[is.na(se.real)]=0
   fixed=rep("",dim(design)[1])
   fixed[fixedparms]="Fixed"
   fixed[boundaryparms]="Boundary"
#
#  Now expand unique parameters to a dataframe with indices
#
   if(replicate.values)
   {
      estimates=NULL
      for (i in 1:dim(data)[1])
      {
         model.indices=model$simplify$pim.translation[indices]
         lookup.indices=match(model.indices,c(rep(0,(i-1)*length(unique(model.indices))),unique(model.indices)))
         covdata=data[rep(i,length(lookup.indices)),]
         xdata=data.frame(vcv.index=lookup.indices,model.index=model.indices,par.index=indices,covdata,estimate=real[lookup.indices],se=se.real[lookup.indices],
               lcl=real.lcl[lookup.indices],ucl=real.ucl[lookup.indices],fixed=fixed[lookup.indices])
         estimates=rbind(estimates,xdata)
      }
      row.names(estimates)=1:dim(estimates)[1]
#
#     Expand v-c matrix
#
      nreals=dim(estimates)[1]
      zz=matrix(1,nrow=nreals,ncol=nreals)*estimates$vcv.index
      vcv.real=matrix(vcv.real[cbind(as.vector(zz),as.vector(t(zz)))],nrow=nreals,ncol=nreals)
   }
   else
   {
      model.indices=model$simplify$pim.translation[indices]
      estimates=data.frame(vcv.index=model.indices,model.index=model.indices,par.index=indices,data,estimate=real,se=se.real,
            lcl=real.lcl,ucl=real.ucl,fixed=fixed)
   }
   reals[[j]]=subset(estimates,select=c("estimate","se"))
#
#  If this is a single model return results
#
   if(number.of.models==1)
   {
      return(list(estimates=estimates,vcv=vcv.real))
   }
   else
#  otherwise construct list for model averaging
   {
	   if(firstmodel)
	   {
		   firstmodel=FALSE
		   nreals=dim(estimates)[1]
		   estimate.mat=matrix(NA,nrow=number.of.models,ncol=nreals)
		   se.mat=matrix(NA,nrow=number.of.models,ncol=nreals)
		   estimates.average=estimates
		   weight=vector("numeric",length=number.of.models)
		   mavg.vcv=vector("list",length=number.of.models)
	   }
	   mavg.vcv[[j]]=vcv.real
	   estimate.mat[j,]=reals[[j]]$estimate
	   if(any(reals[[j]]$se<0)) 
	   {
		   warning("Negative variances for parameters ",paste((1:nreals)[reals[[j]]$se<0],collapse=", ")," for model ",j,". Setting those variances to 0")
		   reals[[j]]$se[reals[[j]]$se<0]=0
	   }
	   se.mat[j,]=reals[[j]]$se
	   weight[j]=model.table$weight[as.numeric(row.names(model.table))==j]	   
   }
#
#  Otherwise if this is the first model in the list setup dataframe for
#  average estimates and average correlation matrix.
#
#   if(firstmodel)
#   {
#      firstmodel=FALSE
#      nreals=dim(estimates)[1]
#      estimates.average=estimates
#      estimates.average$estimate=0
#      cor.average=matrix(0,nrow=nreals,ncol=nreals)
#   }
#
#  For each model the estimates are averaged and so is the correlation matrix
#
#   estimates.average$estimate=estimates.average$estimate+
#         estimates$estimate*model.table$weight[as.numeric(row.names(model.table))==j]
#   cor=vcv.real/outer(estimates$se,estimates$se,"*")
#   if(any(is.infinite(diag(cor)))) 
#	   warning("Infinite correlation (se=0) for model  ",j, " for estimate ",which(is.infinite(diag(cor))),"\n")
#   diag(cor)=1
#   cor.average=cor.average+cor*model.table$weight[as.numeric(row.names(model.table))==j]
}
if(nreals==1){
	estimate.mat=as.vector(estimate.mat)
	se.mat=as.vector(se.mat)
}
if(!mata)
	mavg.res=model.average.list(x=list(estimate=estimate.mat,weight=weight,vcv=mavg.vcv),revised=revised, mata=mata, normal.lm=normal.lm, residual.dfs=residual.dfs, alpha=alpha,...)
else
	mavg.res=model.average.list(x=list(estimate=estimate.mat,weight=weight,se=se.mat),revised=revised, mata=mata, normal.lm=normal.lm, residual.dfs=residual.dfs, alpha=alpha,...)
#
#  After processing each model, the model averaged se and v-c matrix is
#  computed.
#
#revised=TRUE
#se.average=rep(0,nreals)
#for (i in 1:dim(model.table)[1])
#{
#   if(i %in% dropped.models) next   
#   if(revised)
#	   se.average=se.average+model.table$weight[as.numeric(row.names(model.table))==i]*
#			   (reals[[i]]$se^2 + (reals[[i]]$estimate-estimates.average$estimate)^2)  
#   else
#	   se.average=se.average+model.table$weight[as.numeric(row.names(model.table))==i]*
#			   sqrt(reals[[i]]$se^2 + (reals[[i]]$estimate-estimates.average$estimate)^2)
#}
#if(revised) se.average=sqrt(se.average)

names(reals)=1:number.of.models
estimates.average$estimate=mavg.res$estimate
estimates.average$se=mavg.res$se
vcv.real=mavg.res$vcv
#vcv.real=cor.average*outer(se.average,se.average,"*")
#vcv.real[is.nan(vcv.real)]=0
#vcv.real[is.infinite(abs(vcv.real))]=0
if(!mata)
{
	link.list=compute.links.from.reals(estimates.average$estimate,model.list[[1]],parm.indices=estimates.average$par.index,vcv.real=vcv.real,use.mlogits=FALSE)
	estimates.average$lcl=link.list$estimates-1.96*sqrt(diag(link.list$vcv))
	estimates.average$ucl=link.list$estimates+1.96*sqrt(diag(link.list$vcv))
	estimates.average$lcl=apply(data.frame(x=estimates.average$lcl,links=link.list$links),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
	estimates.average$ucl=apply(data.frame(x=estimates.average$ucl,links=link.list$links),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
	estimates.average$lcl[is.na(estimates.average$lcl)]=estimates.average$estimate[is.na(estimates.average$lcl)]
	estimates.average$ucl[is.na(estimates.average$ucl)]=estimates.average$estimate[is.na(estimates.average$ucl)]
}
else
{
	estimates.average$lcl=mavg.res$lcl
	estimates.average$ucl=mavg.res$ucl
}
return(list(estimates=estimates.average,vcv=vcv.real,reals=reals))
}
