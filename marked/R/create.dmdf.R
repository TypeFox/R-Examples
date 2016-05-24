#' Creates a dataframe with all the design data for a particular parameter in a
#' crm model
#' 
#' Creates a dataframe with all the design data for a particular parameter in a
#' crm model which currently includes "cjs" or "js".  These design data are
#' fundamentally different than the design data created for \code{mark} models
#' as explained below.
#' 
#' This function is intended to be called from \code{\link{make.design.data}}.
#' It takes the data in \code{x} and creates a dataframe with all of the data
#' needed to fit capture-recapture models (crm) which currently includes "cjs"
#' (Cormack-Jolly-Seber) or "js" (POPAN formulation of the Jolly-Seber model).
#' Before jumping into the details it is useful to have an understanding of the
#' differences between MARK (via the \code{mark in RMark} function) and the
#' package \code{mra} written by Trent McDonald and how they relate to the
#' implementation in \code{\link{cjs}}.  With MARK, animals can be placed in
#' groups and the parameters for the model specified via PIMs (parameter
#' index matrices) link the parameters to the specific animals.  For
#' example, if for a particular group the \code{Phi} PIM is
#' 
#' \preformatted{ 
#'   1 2 3 
#'     4 5 
#'       6 
#' }
#' 
#' Then animals in that group that were first caught/released on occasion 1
#' have the parameters 1,2,3 for the 3 occasions.  Those first caught/released
#' on occasion 2 have parameters 4 and 5 for occasions 2 and 3 and those first
#' caught/released on occasion 3 have parameter 6 for occasion 3.  Another
#' group of animals would have a different set of indices for the same set of
#' parameters so they could be discriminated.  Many people find this rather
#' confusing until they get used to it but even then if you have many different
#' groups and many occasions, then the indexing process is prone to error.
#' Thus, the rationale for RMark which automates the PIM construction and its
#' use is largely transparent to the user. What RMark does is to create a data
#' structure called design data that automatically assigns \code{design data}
#' to the indices in the PIM.  For example, 1 to 6 would be given the data used
#' to create that group and 1 to 3 would be assigned to cohort 1, ... and 1
#' would be assigned to occasion 1 and 2 and 4 would be assigned to occasion 2
#' etc.  It also creates an age field which follows the diagonals and can be
#' initialized with the intial age at the time of first capture which is group
#' specific.  With a formula and these design data a design matrix can be
#' constructed for the model where the row in the design matrix is the
#' parameter index (e.g., the first 6 rows would be for parameters 1 to 6 as
#' shown above).  That would be all good except for individual covariates which
#' are not group-specific.  MARK handles individual covariates by specifying
#' the covariate name (eg "weight") as a string in the design matrix.  Then for
#' each capture history in plugs in the actual covariate values for that animal
#' to complete the design matrix for that animal. For more details see Laake
#' and Rexstad (2008).
#' 
#' From a brief look at package \code{mra} and personal communication with
#' Trent McDonald, I give the following brief and possibly incorrect
#' description of the pacakge \code{mra} at the time of writing (28 Aug 2008).
#' In that package, the whole concept of PIMS is abandoned and instead
#' covariates are constructed for each occasion for each animal.  Thus, each
#' animal is effectively put in its own group and it has a parameter for each
#' occasion.  This novel approach is quite effective at blurring the lines
#' between design data and individual covariate data and it removes the needs
#' for PIMS because each animal (or unique capture history) has a real
#' parameter for each occasion.  The downside of the pacakge \code{mra} is that
#' every covariate is assumed to be time-varying and any factor variables like
#' \code{time} are coded manually as dummy variables for each level rather than
#' using the R facilities for handling factor variables in the formula to
#' create the design matrix.
#' 
#' In the \code{crm,cjs,js} functions in this package I have used basic idea in
#' \code{mra} but I have taken a different approach to model development that
#' allows for time-varying covariates but does not restrict each covariate to
#' be time-varying and factor variables are used as such which removes the need
#' to construct dummy variables; although the latter could still be used. First
#' an example is in order to explain how this all works. Consider the follwing
#' set of capture histories for small capture-recapture data set with 4 capture
#' occasions:
#' 
#' \preformatted{ 1001 0111 0011 }
#' 
#' To relate the current structure to the concept of PIMS I define the
#' following matrix
#' 
#' \preformatted{
#' 
#'  1 2 3 
#'  4 5 6 
#'  7 8 9 
#' }
#' 
#' If you think of these as \code{Phi} parameter indices, then 1 to 3 are
#' survivals for the intervals 1-2,2-3,3-4 for animal 1, and 4-6 and 7-9 are
#' the same for animals 2 and 3.  This matrix would have a row for each animal.
#' Now you'll notice that for animal 2 parameter 4 is not needed and for animal
#' 3, parameters 7 and 8 are not needed because they are prior to their entry
#' in the study.  While that is certainly true there is no harm in having them
#' and the advantage comes in being able to have a complete matrix in R rather
#' than having a triangular one.
#' 
#' So now we are finally getting to the description of what this function does.
#' It constructs a dataframe with a row for each animal-occasion.  Following on
#' with the example above, depending on how the arguments are set the following
#' dataframe could be constructed:
#' 
#' \preformatted{ row time Time cohort Cohort age Age initial.age 
#'                 1    1    0     1     0     0     0    0
#'                 2    2    1     1     0     1     1    0 
#'                 3    3    2     1     0     2     2    0 
#'                 4    1    0     2     1     0     0    0 
#'                 5    2    1     2     1     1     1    0 
#'                 6    3    2     2     1     2     2    0 
#'                 7    1    0     3     2     0     0    0
#'                 8    2    1     3     2     1     1    0 
#'                 9    3    2     3     2     2     2    0 
#' }
#' 
#' The fields starting with a lowercase character (time,cohort,age) are created
#' as factor variables and those with an uppercase are created as numeric
#' variables.  Note: the \code{age} field is bounded below by the minimum
#' \code{initial.age} to avoid creating factor levels with non-existent data
#' for occasions prior to first capture that are not used. For example, an
#' animal first caught on occasion 2 with an \code{initial.age=0} is
#' technically -1 on occasion 1 with a \code{time.interval} of 1. However, that
#' parameter would never be used in the model and we don't want a factor level
#' of -1.
#' 
#' A formula of ~time would create a design matrix with 3 columns (one for each
#' factor level) and ~Time would create one with 2 columns with the first being
#' an intercept and the second with the numeric value of Time.
#' 
#' Now here is the simplicity of it.  The following few expressions in R will
#' convert this dataframe into a matrix of real parameters (assuming
#' \code{beta=c(1,1,1)} that are structured like the square PIM matrix without
#' the use of PIMs.
#' 
#' \preformatted{ 
#' nocc=4
#' x=data.frame(ch=c("1001","0111","0011"),stringsAsFactors=FALSE)
#' beta=c(1,1,1) 
#' x.proc=process.data(x,model="cjs")
#' Phi.dmdf=make.design.data(x.proc)$Phi 
#' Phi.dm=create.dm(Phi.dmdf,~time)
#' Phimat=matrix(plogis(Phi.dm%*%beta),nrow=nrow(x),ncol=nocc-1,byrow=TRUE) 
#' }
#' 
#' Note that the order of the columns for \code{Phi.dmdf} differs slightly from
#' what is shown above. Also, \code{plogis} is an R function that computes the
#' inverse-logit. Once you have the matrix of \code{Phi} and \code{p} values
#' the calculation of the likelihood is straightforward using the formulation
#' of Pledger et al. (2003) (see \code{\link{cjs.lnl}}). The values in the
#' design dataframe are not limited to these fields.  The 2 arguments
#' \code{time.varying} and \code{fields} are vectors of character strings which
#' specify the names of the dataframe columns in \code{x} that should be
#' included. For example if \code{x} contained a field \code{sex} with the
#' values "M","F","M" for the 3 records in our example, and the argument
#' \code{fields=c("sex")} was used then a column named \code{sex} would be
#' included in design dataframe with the values
#' "M","M","M","F","F","F","M","M","M".  The value of the column \code{sex} in
#' \code{x} is repeated for each of the occasions for that animal(capture
#' history). Now if the value of the field changes for each occasion then we
#' use the argument \code{time.varying} instead.  To differentiate the values
#' in the dataframe \code{x} the columns are named with an occasion number.
#' For example, if the variable was named \code{cov} and it was to be used for
#' \code{Phi}, then the variables would be named \code{cov1,cov2,cov3} in
#' \code{x}. Let's say that x was structured as follows:
#' 
#' \preformatted{ 
#' ch   cov1 cov2 cov3 
#' 1001   1   0     1 
#' 0111   0   2     1 
#' 0011   0   0     0 
#' }
#' 
#' If you specified the argument \code{time.varying=c("cov")} then in the
#' design dataframe a field named \code{cov} would be created and the values
#' would be \code{1,0,1,0,2,1,0,0,0}. Thus the value is both animal and
#' occasion specific.  Had the covariate been used for \code{p} then they would
#' be named \code{cov2,cov3,cov4} because the covariate is for those occasions
#' for \code{p} whereas for \code{Phi} the covariate is labelled with the
#' occasion that begins the interval.  Any number of fields can be specified in \code{fields} and
#' \code{time.varying} that are specified in \code{x}.
#' 
#' The input dataframe \code{x} has a few minor requirements on its structure.
#' First, it must contain a field called \code{ch} which contains the
#' capture-history as a string.  Note that in general strings are converted to
#' factor variables by default when they are put into a dataframe but as shown
#' above that can be controlled by the argument \code{stringsAsFactors=FALSE}.
#' The capture history should be composed only of 0 or 1 and they should all be
#' the same length (at present no error checking on this). Although it is not
#' necessary, the dataframe can contain a field named \code{freq} which
#' specifies the frequency of that capture history.  If the value of
#' \code{freq} is negative then these are treated as loss on capture at the
#' final capture occasion (last 1).  If \code{freq} is missing then a value of
#' 1 is assumed for all records.  Another optional field is \code{initial.age}
#' which specifies the age of the animal at the time it was first captured.
#' This field is used to construct the \code{age} and \code{Age} fields in the
#' design dataframe.  The default is to assume \code{initial.age=0} which means
#' the \code{age} is really time since first marked.  Any other fields in
#' \code{x} are user-specified and can be a combination of factor and numeric
#' variables that are either time-invariant or time-varying (and named
#' appropriately).
#' 
#' The behavior of create.dmdf can vary depending on the values of begin.time and 
#' time.intervals. An explanation of these values and how they can be used is given 
#' in \code{\link{process.data}}.
#' 
#' @param x processed dataframe from function \code{\link{process.data}}
#' @param parameter list with fields defining each values for each parameter;
#' as created by \code{\link{setup.parameters}}
#' @param time.varying vector of field names that are time-varying for this
#' parameter
#' @param fields character vector containing field names for variables in x to
#' be included in design matrix dataframe; if NULL all other than ch are
#' included
#' @return A dataframe with all of the individual covariate data and the
#' standard design data of time, Time, cohort, Cohort, age and Age; where
#' lowercase first letter implies a factor and uppercase is a numeric variable
#' for those variables.
#' @author Jeff Laake
#' @references Laake, J. and E. Rexstad (2007). RMark -- an alternative
#' approach to building linear models in MARK. Program MARK: A Gentle
#' Introduction. E. Cooch and G. C. White.
#' 
#' Pledger, S., K. H. Pollock, et al. (2003). Open capture-recapture models
#' with heterogeneity: I. Cormack-Jolly-Seber model. Biometrics 59(4):786-794.
create.dmdf=function(x,parameter,time.varying=NULL,fields=NULL)
{
   last= -parameter$num
   begin.num=parameter$begin+1
   chp = process.ch(x$data$ch)
   firstseen=chp$first
   nocc=x$nocc
   time.intervals=x$time.intervals
   if(is.null(x$data$begin.time))
	   begin.time=x$begin.time
   else
	   begin.time=x$data$begin.time
#  create base df
   df=create.base.dmdf(x,parameter)
   if(is.null(df))return(NULL)
#  add time,age,choort
   if(length(begin.time)==1 & is.vector(time.intervals))
   {
	   timedf=data.frame(occ=1:nocc,time=begin.time+c(0,cumsum(time.intervals)))
	   if(parameter$interval)
		   if(!all(time.intervals==1))
		      timedf$time.interval=c(time.intervals,NA)
	   df=merge(df,timedf,by="occ")
	   timedf=data.frame(id=x$data$id,cohort=timedf$time[firstseen])
	   df=merge(df,timedf,by="id")
   }else
   {
	   if(length(begin.time)==1) 
		   begin.time=rep(begin.time,nrow(x$data))
	   if(is.vector(time.intervals)) 
		   time.intervals=matrix(time.intervals,nrow=nrow(x$data),ncol=length(time.intervals),byrow=TRUE)
	   cum.time.intervals=t(apply(time.intervals,1,cumsum))
	   times=begin.time+cbind(rep(0,nrow(x$data)),cum.time.intervals)
	   cohort=rep(times[cbind(1:nrow(x$data),firstseen)],each=ncol(times))
	   if(!parameter$interval | all(time.intervals==1))
	       timedf=data.frame(occid=apply(expand.grid(occ=1:nocc,id=factor(1:nrow(x$data))),1,paste,collapse=""),
			       cohort=cohort,time=as.vector(t(times)))
	   else
		   timedf=data.frame(occid=apply(expand.grid(occ=1:nocc,id=factor(1:nrow(x$data))),1,paste,collapse=""),
				   cohort=cohort,time=as.vector(t(times)),time.interval=as.vector(t(cbind(time.intervals,rep(NA,nrow(x$data))))))
	   df$occid=paste(df$occ,df$id,sep="")
	   df=merge(df,timedf,by="occid")
	   df$idocc=NULL
   }
   df=df[order(df$seq),]
   df$age=df$time-df$cohort
   if(!is.null(x$data$initial.age)) df$age=df$age+x$data$initial.age[df$id]
#  if there are any time varying covariates then construct the data matrix for
#  those covariates. That matrix is appended to other data in x that is non-time varying
   times=df$time[df$id==1]
 #  occ=1:(nocc+parameter$num)
   tcv=NULL
   if(!is.null(time.varying))
   {
	   for (i in 1:length(time.varying))
	   {
		   vnames=paste(time.varying[i],times,sep="")     
		   if( !all(vnames %in% names(x$data)))
			   stop("Missing time varying variable ",paste(vnames[!vnames%in%names(x$data)],collapse=","))
		   if(i==1) 
			   tcv=data.frame(as.vector(t(as.matrix(x$data[, vnames]))))
		   else
			   tcv=cbind(tcv,as.vector(t(as.matrix(x$data[, vnames]))))
		   x$data=x$data[,!names(x$data)%in%vnames]
	   }  
	   names(tcv)=time.varying
   }
   if(!is.null(tcv))  
      df=cbind(df,tcv)
 # add static fields
   if(is.null(fields))
	   fields=names(x$data)[!names(x$data)%in%c("ch","initial.age")]
   else
	   fields=c(fields,"id")
   df=merge(df,x$data[,fields],by="id") 
   df=df[order(df$seq),]
   df$seq=NULL
   df$Time=df$time-min(df$time)
   df$Cohort=df$cohort-min(df$cohort)
  # df$Age=df$age-max(0,min(df$age))
   df$age[df$age<0]=0
   df$time=factor(df$time)
   df$Age=df$age
   df$age=factor(df$age) 
   df$cohort=factor(df$cohort)
   if("group"%in%names(df))
	  levels(df$group)=apply(x$group.covariates,1,paste,collapse="")
   if(!is.null(x$strata_data)&!is.null(df$stratum))
	   df=cbind(df,x$strata_data)
   return(df)
}

create.base.dmdf=function(x,parameter)
{
	last= -parameter$num
	begin.num=parameter$begin+1
	nocc=x$nocc + parameter$num
	occasions=begin.num:(parameter$begin+nocc)
	sl=factor(x$strata.labels)
	nstrata=length(sl)
    # MVMS model is currently split off but eventually this should be merged in with
	# other bi-level MS models
	if(nchar(x$model)>=4 & substr(x$model,1,4)=="MVMS")
	{
		if(!is.null(parameter$obs) & parameter$obs)
			if(any(x$strata.list$uncertain))
			   dfl=mvms_design_data(x$strata.list$df.states,x$strata.list$df,transition=parameter$tostrata)
	        else
				return(NULL)
		else
		   dfl=mvms_design_data(x$strata.list$df.states,transition=parameter$tostrata)
		df=expand.grid(occ=occasions,id=factor(1:nrow(x$data)))
		dfl=dfl[rep(1:nrow(dfl),each=nrow(df)),,drop=FALSE]
		df=cbind(df,dfl)
		if(parameter$tostrata)
			df=df[order(df$id,df$occ,df$stratum,df$tostratum),]
		else
			df=df[order(df$id,df$occ,df$stratum),]
		df$seq=1:nrow(df)	
		rownames(df)=df$seq
	} else
	{
		if(!is.null(x$strata.list))
		{
			st.index=which("states"==names(x$strata.list))
			oth.index=which("states"!=names(x$strata.list))
			oth.name=names(x$strata.list)[oth.index]
			oth=factor(rep(x$strata.list[[oth.index]],each=length(x$strata.list$states)))
			states=factor(rep(x$strata.list$states,times=length(x$strata.list[[oth.index]])))
			# For 2-level independent strata model 
			if(parameter$whichlevel!=0)
			{
				if(parameter$whichlevel==1)
					sl=factor(x$strata.list$states)
				else
					sl=factor(x$strata.list[[oth.index]])
				nstrata=length(sl)
				x$strata.list=NULL
			}
		}
		# not by stratum
		if(!parameter$bystratum)
		{
			df=expand.grid(occ=occasions,id=factor(1:nrow(x$data)))
			df$seq=1:nrow(df)
		}
		else
		{
			# from stratum - to stratum transition
			if(parameter$tostrata)
			{
				df=expand.grid(tostratum=sl[1:nstrata],stratum=sl[1:nstrata],occ=occasions,id=factor(1:nrow(x$data)))
				df=df[,c("stratum","tostratum","occ","id")]
				df$seq=1:nrow(df)
			}
			# by stratum
			else
			{
				df=expand.grid(stratum=sl[1:nstrata],occ=occasions,id=factor(1:nrow(x$data)))
				df$seq=1:nrow(df)
			}
			# two-level stratification
			if(!is.null(x$strata.list))
			{
				state.df=data.frame(stratum=sl[1:nstrata],state=states,oth=oth)	
				df=merge(df,state.df,by="stratum")
				if(parameter$tostrata)
				{
					state.df=data.frame(tostratum=sl[1:nstrata],tostate=states,tooth=oth)	
					df=merge(df,state.df,by="tostratum")
				} 
				if(parameter$whichlevel!=0)
				{
					if(parameter$whichlevel==1){
						names(df)[names(df)=="stratum"]="state"
						names(df)[names(df)=="tostratum"]="tostate"
					} 
					if(parameter$whichlevel==2) {
						names(df)[names(df)=="stratum"]=oth.name
						names(df)[names(df)=="tostratum"]=paste("to",oth.name,sep="")
					}	 
				}else
				{
					names(df)[names(df)=="oth"]=oth.name
					names(df)[names(df)=="tooth"]=paste("to",oth.name,sep="")
				}
			}else
			{
				if(parameter$whichlevel!=0)
				{
					if(parameter$whichlevel==1)names(df)[names(df)=="stratum"]="state"
					if(parameter$whichlevel==2) names(df)[names(df)=="stratum"]=oth.name
					if(parameter$whichlevel==1)names(df)[names(df)=="tostratum"]="tostate"
					if(parameter$whichlevel==2) names(df)[names(df)=="tostratum"]=paste("to",oth.name,sep="")
				}	 
			}
		}
	}
	return(df)
}




