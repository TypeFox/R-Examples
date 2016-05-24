#' Compute design data for a specific parameter in the MARK model (internal
#' use)
#' 
#' For a specific type of parameter (e.g., Phi, p, r etc), it creates a data
#' frame containing design data for each parameter of that type in the model as
#' structured by an all different PIM (parameter information matrix). The
#' design data are used in constructing the design matrix for MARK with
#' user-specified model formulae as in \code{\link{make.mark.model}}.
#' 
#' This function is called by \code{\link{make.design.data}} to create all of
#' the default design data for a particular type of model and by
#' \code{\link{add.design.data}} to add binned design data fields for a
#' particular type of parameter. The design data created by this function
#' include \code{group}, \code{age}, \code{time} and \code{cohort} as factors
#' variables and continuous (non-factor) versions of all but \code{group}.  In
#' addition, if groups have been defined for the data, then a data column is
#' added for each factor variable used to define the groups.  Also for specific
#' closed capture heterogeneity models (\code{model}="HetClosed", "FullHet",
#' "HetHug", "FullHetHug") the data column \code{mixture} is added to the
#' design data. The arguments for this function are defined for each model by
#' the function \code{\link{setup.model}}.
#' 
#' @param data data list created by \code{\link{process.data}}
#' @param begin 0 for survival type, 1 for capture type
#' @param num number of parameters relative to number of occasions (0 or -1)
#' @param type type of parameter structure (Triang (STriang) or Square)
#' @param mix if TRUE this is a mixed parameter
#' @param rows number of rows relative to number of mixtures
#' @param pim.type type of pim structure; either all (all-different) or time
#' @param secondary TRUE if a parameter for the secondary periods of robust
#' design
#' @param nstrata number of strata for multistrata
#' @param tostrata set to TRUE for transition parameters
#' @param strata.labels labels for strata as identified in capture history
#' @param subtract.stratum for each stratum, the to.strata that is computed by
#' subtraction
#' @param common.zero if TRUE, uses a common begin.time to set origin (0) for
#' Time variable defaults to FALSE for legacy reasons but should be set to TRUE
#' for models that share formula like p and c with the Time model
#' @param sub.stratum the number of strata to subtract for parameters that use 
#' mlogit across strata like pi and Omega for RDMSOpenMisClass
#' @param limits For RDMSOccRepro values that set row and col (if any) start on states
#' @return design.data: a data frame containing all of the design data fields
#' for a particular type of parameter \item{group}{group factor level}
#' \item{age}{age factor level} \item{time}{time factor level}
#' \item{cohort}{cohort factor level} \item{Age}{age as a continuous variable}
#' \item{Time}{time as a continuous variable} \item{Cohort}{cohort as a
#' continuous variable} \item{mixture}{mixture factor level} \item{other
#' fields}{any factor variables used to define groups}
#' @author Jeff Laake
#' @seealso \code{\link{make.design.data}}, \code{\link{add.design.data}}
#' @keywords utility
compute.design.data <-
function(data,begin,num,type="Triang",mix=FALSE,rows=0,pim.type="all",
           secondary,nstrata=1,tostrata=FALSE,strata.labels=NULL,
           subtract.stratum=strata.labels,common.zero=FALSE,sub.stratum=0,limits=NULL)
{
# -------------------------------------------------------------------------------------------------------------
#
# compute.design.data -  creates a design dataframe that is used to construct the design matrix 
#
# Value:
#
#  design.data     - design data (cohort, group, age, time) for a particular parameter
#
#
# -------------------------------------------------------------------------------------------------------------
#
# Create a data matrix for the parameter PIM that contains age, year and cohort for each index
# This data matrix (design.data) is used to create the design matrix from the formulas
#
  if(secondary)
     num.sessions=data$nocc
  else
  {
     num.sessions=1
     if(!is.null(pim.type)&&pim.type=="constant")
       num=1
     else
       num=data$nocc+num
  }
  if(type=="Square")
  {
     num.lines=1
     if(is.null(mix) || !mix)
        num.rows=1
     else
        num.rows=data$mixtures+rows 
  }
  else
  {
  #
  #  pim.type field allows either all-different or time pims for Triangular pims
  #
	 if(is.null(mix) || !mix)
		num.rows=1
	 else
		num.rows=data$mixtures+rows 
     if(pim.type=="all")
     {
        num.lines=num
     }
     else
     {
        num.lines=1
     }
  }
  if(setup.model(data$model,data$nocc)$robust)
     time.intervals=data$time.intervals[data$time.intervals>0]
  else
     time.intervals=data$time.intervals
  number.of.groups=dim(data$freq)[2]
  design.data=NULL
  nsubtract.stratum=match(subtract.stratum,strata.labels)
  all.tostrata=FALSE
  if(sub.stratum==-1)
  {
	  all.tostrata=TRUE
	  sub.stratum=0
  }
  if(is.null(limits) || limits[1]=="0")
	  start.stratum=1
  else
	  start.stratum=as.numeric(limits[1])
  for(j in 1:number.of.groups)
  for (jj in start.stratum:(nstrata-sub.stratum))
  for(l in 1:num.sessions)
  {
      if(tostrata)
	  {
		 if(!all.tostrata)
		 {
			 if(is.null(limits))
				 other.strata= sequence(nstrata)[sequence(nstrata)!=nsubtract.stratum[jj]]
			 else
			 {
				 if(limits[1]=="0")
					 other.strata=(as.numeric(limits[2])+1):nstrata
				 else
					 other.strata=jj:nstrata
			 }		 
		 }
	     else
			 other.strata= 1:nstrata
	  }		 
      else
         other.strata=1
  for(to.strata in other.strata)
  {
     if(secondary)
     {
       if(is.na(num))
          ncol=1
       else
       {
          ncol=data$nocc.secondary[l]+num
#		  if(any(ncol<1))ncol=rep(1,length(data$nocc.secondary))
          if(type%in%c("Triang","STriang"))num.lines=ncol
       }
     }
     else
        ncol=num
	 ncol.save=ncol
	 for(k in 1:num.rows)
	 {
	    ncol=ncol.save
  	    for(i in 1:num.lines)
        {    
#
#        Define age variable
#
        if(secondary)
           ages=0
        else
           if(begin==0)
              ages=data$initial.ages[j]+data$age.unit*(cumsum(c(0,time.intervals[i:num]))[1:ncol])
           else
              ages=data$initial.ages[j]+data$age.unit*(cumsum(time.intervals[i:num]))
#
#     Define cohort variable
#
        if(secondary)
           if(!type%in%c("Triang","STriang"))
              cohort=0
           else
              cohort=i
        else
          if(i==1)
            if(length(data$begin.time)==1)
               cohort=data$begin.time
            else
               cohort=data$begin.time[j]
         else
            if(length(data$begin.time)==1)
               cohort=data$begin.time+sum(time.intervals[1:(i-1)])
            else
               cohort=data$begin.time[j]+sum(time.intervals[1:(i-1)])
#
#     Define time variable
#
        if(secondary)
          if(is.na(num))
             times=0
          else
             if(type%in%c("Triang","STriang"))
                times=(begin+i):(data$nocc.secondary[l]+num)
             else
                times=(begin+1):(begin+ncol)
        else
          if(begin==0)
             if(i==num)
                times=cohort
             else
                times=c(cohort,cohort+cumsum(time.intervals[i:(num-1)]))
          else
             times=cohort+cumsum(time.intervals[i:num])
#
#      Create design data as needed for the parameter
#
        if(type%in%c("Triang","STriang"))
        {
          if(pim.type=="all")
          {
             add.design.data=cbind(rep(j,ncol),rep(cohort,ncol),ages,times,(i-1)+1:ncol,rep(i,ncol))
             dd.names=c("group","cohort","age","time","occ","occ.cohort")
          }
          else
            if(pim.type=="time")
            {
                add.design.data=cbind(rep(j,ncol),times)
                dd.names=c("group","time")
            }
            else
				if(pim.type=="age")
				{
					add.design.data=cbind(rep(j,ncol),rep(cohort,ncol),ages,times)
					dd.names=c("group","cohort","age","time")
				}
				else	
				{
                add.design.data=matrix(rep(j,ncol),nrow=1)
                dd.names=c("group")
                }
        }
        else
        {
          add.design.data=cbind(rep(j,ncol),ages,times)
          dd.names=c("group","age","time")
        }
        if(!is.null(mix) && mix)
        {
          add.design.data=cbind(add.design.data,rep(k,ncol))
          dd.names=c(dd.names,"mixture")
        }
        if(nstrata>1|tostrata)
           if(tostrata)
           {
              add.design.data=cbind(add.design.data,rep(jj,ncol),rep(to.strata,ncol))
              dd.names=c(dd.names,"stratum","tostratum")
           }
           else
           {
              add.design.data=cbind(add.design.data,rep(jj,ncol))
              dd.names=c(dd.names,"stratum")
           }
        if(secondary)
        {
          add.design.data=cbind(add.design.data,rep(l,ncol))
          dd.names=c(dd.names,"session")
        }
#
#     Add rows to existing design data
#
#        if(is.null(design.data) || ncol(add.design.data)==ncol(design.data))
          design.data=rbind(design.data,add.design.data)
#		else
#		  browser()
#
#      If trianular pim type, decrement number of cols
#
       if(type%in%c("Triang","STriang"))
          ncol=ncol-1
       }
	 }
  }
  }
   design.data=as.data.frame(design.data,row.names=NULL)
   names(design.data)=dd.names
#
#  Add Cohort, Age and Time variables
#
   if(!is.null(design.data$cohort))design.data$Cohort=design.data$cohort- min(design.data$cohort)
   if(!is.null(design.data$age))design.data$Age=design.data$age
   if(!is.null(design.data$time))
     if(common.zero)
        design.data$Time=design.data$time- data$begin.time     
     else
        design.data$Time=design.data$time- min(design.data$time)
   if(nstrata>1)
      design.data$stratum=as.factor(strata.labels[design.data$stratum])
      if(!is.null(design.data$tostratum))
        design.data$tostratum=as.factor(strata.labels[design.data$tostratum])
#
#  Next add grouping variables
#
   if(!is.null(data$group.covariates))
   {
      ix=grep("age",names(data$group.covariates))
      cnames=names(data$group.covariates)
      if(length(ix)!=0)
         if(names(data$group.covariates)[ix]=="age")
         {
            cnames[ix]="initial.age.class"
            names(data$group.covariates)=cnames
         }
      gc=data.frame(data$group.covariates[design.data$group,]) 
      names(gc)=cnames 
      row.names(gc)=NULL
      design.data=cbind(design.data,gc)
   }
#
#  Finally if there are stratum variables, add dummy variables for each
#
   if(!is.null(design.data$stratum))
      for (label in strata.labels)
      {
         design.data[label]=0
         design.data[design.data$stratum==label,label]=1
      }
   if(!is.null(design.data$tostratum))
      for (label in strata.labels)
      {
         design.data[paste("to",label,sep="")]=0
         design.data[design.data$tostratum==label,paste("to",label,sep="")]=1
      }
#  Remove occ field unless this is a Multistrata model
   if(data$model!="Multistrata")design.data$occ=NULL
   return(design.data)
}
