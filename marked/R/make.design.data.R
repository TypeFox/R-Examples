#' Create design dataframes for crm
#' 
#' For each type of parameter in the analysis model (e.g, p, Phi), this
#' function creates design data for model fitting from the original data. 
#' These design data expand the original data record for a single(freq=1)/multiple(freq>1) individuals
#' to many records where each record is for a single occasion. The design data can
#' be specific to the parameter and a list of design data frames are returned with a dataframe
#' for each parameter.
#' 
#' For each record in the design data, default data fields are created based on the model. For example,
#' for CJS, the default fields include cohort, age, time which are factor variables and Cohort, Age, and Time
#' which are numeric versions of the same fields.  The values of these can be altered by values of 
#' begin.time, time.intervals and initial.ages set in the call to \code{\link{process.data}}. In addition if groups are identifed the
#' grouping factor variables are added to the design data. Any other fields in the data are repeated on each record
#' for the animal(s) (eg weight), unless they are defined as time.varying in which case the fields should be named
#' with the convention xn where x is the base name for the variable and n is the time value (eg, td1999, td2000,...).
#' For time varying fields the variable name in the design data is the base name and the value for the record is
#' the one for that occasion(time). The variables can be refined for each parameter by including argument static with the
#' character vector of static fields to include.
#' 
#' After creating design data, you can create a field in the dataframe for a parameter named fix and it can be assigned values
#' at the real parameter for that occasion/id will be fixed at the value. For parameters that are not supposed to be fixed, the field
#' fix should be assigned NA.  For example, ddl$Phi$fix=ifelse(ddl$Phi$time==2,1,NA) will assign all real Phi values to 1 for the interval
#' beginning at time 2. If there is no field fix, then no parameters are fixed.  For mlogit parameters, a fix field is added automatically and
#' the value 1 is specified for the stratum that is supposed to be computed by subtraction and the remainder are set to NA.  If some Psi are not possible
#' then those values can be changed to 0.
#' 
#' The following variable names are reserved and should not be used in the data:
#' occ,age,time,cohort,Age,Time,Cohort,Y,Z,initial.age,begin.time,time.interval,fix
#' 
#' @param data Processed data list; resulting value from process.data
#' @param parameters Optional list containing a list for each type of parameter
#' (list of lists); each parameter list is named with the parameter name (eg
#' Phi); each parameter list can contain vectors named age.bins,time.bins and
#' cohort.bins, time.varying, and static \tabular{ll}{ \code{age.bins} \tab bins
#' for binning ages\cr \code{time.bins} \tab bins for binning times\cr
#' \code{cohort.bins} \tab bins for binning cohorts\cr \code{time.varying} \tab
#' vector of field names that are time varying for this parameter\cr
#' \code{static} \tab vector of field names that are to be included that are
#' not time varying for this parameter\cr }
#' @return The function value is a list of data frames. The list contains a
#' data frame for each type of parameter in the model (e.g., Phi and p for
#' CJS).  The names of the list elements are the parameter names (e.g., Phi).
#' The structure of the dataframe depends on the calling arguments and the
#' model & data structure. In addition the value of parameters argument is saved as design.parameters.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{process.data}},\code{\link{merge_design.covariates}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' dipper.proc=process.data(dipper)
#' ddl=make.design.data(dipper.proc)
#' }
"make.design.data" <-                                                      
function(data,parameters=list())
{
#
# Check validity of parameter list; stop if not valid
#
  if(!valid.parameters(data$model,parameters)) stop()
#
#  Add following elements based on type of model
#           begin            - index for compute.design.data
#           num              - number of parameters relative to number of occasions
#
  par.list=setup.parameters(data$model,check=TRUE)
  parameters=setup.parameters(data$model,parameters,data$nocc,check=FALSE,
          number.of.groups=dim(data$freq)[2])
  parameters=parameters[par.list]
  model.list=setup.model(data$model,data$nocc,data$mixtures)
# Create data for the each parameter in the model with age, year and cohort for each index
# This data matrix (design.data) is used below to create the design matrix from the formulas
# If age,cohort or year bins are given, use those.  Otherwise each is treated as a factor 
# wihtout binning.
full.design.data=vector("list",length=length(parameters))
#
# Create design data for crm models
#  Loop over each parameter in the model
   for(i in 1:length(parameters))
   {
#    Compute design data for parameters other than N
     if(names(parameters)[i]!="N")
	 {
		 if(is.null(parameters[[i]]$time.varying))
			 time.varying=NULL
		 else
			 time.varying=parameters[[i]]$time.varying
		 if(is.null(parameters[[i]]$static))
			 fields=NULL
		 else
			 fields=parameters[[i]]$static   
		 full.design.data[[i]]=create.dmdf(data,parameters[[i]],time.varying=time.varying,fields=fields)
		 if(is.null(full.design.data[[i]]))next
	 }else
#    Compute design data for N
	 {
	    if(!is.null(data$group.covariates))
        {
           xlist=split(data$data,data$data[,names(data$group.covariates),drop=FALSE])
           xlist=xlist[as.vector(sapply(xlist,function(x) nrow(x)))>0]
        } else
        {                              
           xlist=data$data
        }
        numvar=sapply(data$data,is.numeric)
        numvar=numvar[names(numvar)!="freq"]
        if(any(numvar))
        {
           numvar=names(data$data[,names(data$data)!="freq"])[numvar]
           if(!is.null(data$group.covariates))
           {
              xmeans=sapply(xlist,function(x) sapply(subset(x,select=numvar),mean))
			  if(length(numvar)==1)
			  {
				  dd=data.frame(group=1:nrow(data$group.covariates))
				  dd[,numvar]=xmeans
			  }else
			  {
				  dd=data.frame(group=1:nrow(data$group.covariates))
				  dd=cbind(dd,t(xmeans))
			  }
              full.design.data[[i]]=merge(cbind(data.frame(group=1:nrow(data$group.covariates),data$group.covariates)),
					                         dd,by.x="group",all.x=TRUE)
#             full.design.data[[i]]=cbind(data$group.covariates,group=factor(1:nrow(data$group.covariates)),xmeans)
              row.names(full.design.data[[i]])=NULL
              names(full.design.data[[i]])=c("group",names(data$group.covariates),numvar[numvar!="group"])          
           }
           else
           {
              xmeans=colMeans(subset(data$data,select=numvar))         
              if(ncol(t(xmeans))==0)
                 full.design.data[[i]]=data.frame(N=1)
              else
              {
                 full.design.data[[i]]=data.frame(t(xmeans))  
                 row.names(full.design.data[[i]])=NULL
                 names(full.design.data[[i]])=numvar          
              }
           }
         }else
         {
           numvar=NULL
           if(!is.null(data$group.covariates))
           {
              full.design.data[[i]]=data$group.covariates
              full.design.data[[i]]$group=factor(row.names(data$group.covariates))
           }
           else
              full.design.data[[i]]=data.frame(N=1)
         }     
      } 
#	  if(!toupper(data$model)%in%c("PROBITCJS","NULL"))
#		  if("Y" %in% names(full.design.data[[i]]))
#			  full.design.data[[i]]$Y=NULL
	  # assign subtract.stratum and fix values to 1 unless subtract.stratum=="NONE"
      # the code now handles parameters for 2iMSCJS where strata are in levels (eg states,areas)
      labels=data$strata.labels
	  if(parameters[[i]]$whichlevel==1)labels=data$strata.list$states
	  if(parameters[[i]]$whichlevel==2)
	  {
		  oth.index=which("states"!=names(data$strata.list))
		  labels=data$strata.list[[oth.index]]
	  } 
	  if(!is.null(parameters[[i]]$tostrata) && parameters[[i]]$tostrata)
	  {
		  if(is.null(parameters[[i]]$subtract.stratum)) 
		  {
			  if(parameters[[i]]$whichlevel==0)
			  {
				  field="stratum"
			  	  parameters[[i]]$subtract.stratum=data$strata.labels
			  }
			  else
			      if(parameters[[i]]$whichlevel==1)
				  {
					  field="state"  
					  parameters[[i]]$subtract.stratum=data$strata.list$states
				  }
				  else
				  {
					  field=names(data$strata.list)[oth.index]
					  parameters[[i]]$subtract.stratum=data$strata.list[[oth.index]]
				  }
		  }	  	  
		  if(toupper(parameters[[i]]$subtract.stratum)[1]!="NONE")
		  {
			  full.design.data[[i]]$fix=NA
			  if(parameters[[i]]$whichlevel==0&length(parameters[[i]]$subtract.stratum)!=length(data$strata.labels)) stop("\nlength of subtract.stratum does not match number of strata")
			  for(j in 1:length(labels))
			  {
				  if(parameters[[i]]$whichlevel==0&&!parameters[[i]]$subtract.stratum[j]%in%data$strata.labels) stop("\n invalid value of subtract.stratum: ",parameters[[i]]$subtract.stratum[j])
				  full.design.data[[i]]$fix[full.design.data[[i]][field]==labels[j] & full.design.data[[i]][paste("to",field,sep="")]==parameters[[i]]$subtract.stratum[j]]=1
			  } 
		  }		  
	  }
	  # assign p/delta for CJS type models conditionining on first release
      if(parameters[[i]]$cjs)
	  {
		  full.design.data[[i]]$fix=NA
		  full.design.data[[i]]$fix[as.character(full.design.data[[i]]$time)==as.character(full.design.data[[i]]$cohort)]=1	  
	  }	
	  full.design.data[[i]]=droplevels(full.design.data[[i]])
	  if(names(parameters)[i]=="tau"&!is.null(full.design.data[[i]]$tag1))
	  {
		  full.design.data[[i]]$fix=NA
		  full.design.data[[i]]$fix[full.design.data[[i]]$tag1==0&full.design.data[[i]]$tag2==0]=1	  
	  }
	  full.design.data[[i]]$order=1:nrow(full.design.data[[i]])
   }
   names(full.design.data)=names(parameters)
   full.design.data[["design.parameters"]]=parameters
   full.design.data$ehmat=data$ehmat
   return(full.design.data)
}

        
         
