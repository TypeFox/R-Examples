#' Process encounter history dataframe for MARK analysis
#' 
#' Prior to analyzing the data, this function initializes several variables
#' (e.g., number of capture occasions, time intervals) that are often specific
#' to the capture-recapture model being fitted to the data.  It also is used to
#' 1) define groups in the data that represent different levels of one or more
#' factor covariates (e.g., sex), 2) define time intervals between capture
#' occasions (if not 1), and 3) create an age structure for the data, if any.
#' 
#' For examples of \code{data}, see
#' \code{\link{dipper}},\code{\link{edwards.eberhardt}},\code{\link{example.data}}.
#' The structure of the encounter history and the analysis depends on the
#' analysis model to some extent. Thus, it is necessary to process a dataframe
#' with the encounter history (\code{ch}) and a chosen \code{model} to define
#' the relevant values.  For example, number of capture occasions (\code{nocc})
#' is automatically computed based on the length of the encounter history
#' (\code{ch}) in \code{data}; however, this is dependent on the type of
#' analysis model.  For models such as "CJS", "Pradel" and others, it is simply
#' the length of \code{ch}.  Whereas, for "Burnham" and "Barker" models,the
#' encounter history contains both capture and resight/recovery values so
#' \code{nocc} is one-half the length of \code{ch}. Likewise, the number of
#' \code{time.intervals} depends on the model.  For models, such as "CJS",
#' "Pradel" and others, the number of \code{time.intervals} is \code{nocc-1};
#' whereas, for capture&recovery(resight) models the number of
#' \code{time.intervals} is \code{nocc}. The default time interval is unit time
#' (1) and if this is adequate, the function will assign the appropriate
#' length.  A processed data frame can only be analyzed using the model that
#' was specified.  The \code{model} value is used by the functions
#' \code{\link{make.design.data}}, \code{\link{add.design.data}}, and
#' \code{\link{make.mark.model}} to define the model structure as it relates to
#' the data. Thus, if the data are going to be analysed with different
#' underlying models, create different processed data sets with the model name
#' as an extension.  For example, \code{dipper.cjs=process.data(dipper)} and
#' \code{dipper.popan=process.data(dipper,model="POPAN")}.
#' 
#' This function will report inconsistencies in the lengths of the capture
#' history values and when invalid entries are given in the capture history.
#' For example, with the "CJS" model, the capture history should only contain 0
#' and 1 whereas for "Barker" it can contain 0,1,2.  For "Multistrata" models,
#' the code will automatically identify the number of strata and strata labels
#' based on the unique alphabetic codes used in the capture histories.
#' 
#' The argument \code{begin.time} specifies the time for the first capture
#' occasion.  This is used in creating the levels of the time factor variable
#' in the design data and for labelling parameters. If the \code{begin.time}
#' varies by group, enter a vector of times with one for each group. Note that
#' the time values for survivals are based on the beginning of the survival
#' interval and capture probabilities are labeled based on the time of the
#' capture occasion.  Likewise, age labels for survival are the ages at the
#' beginning times of the intervals and for capture probabilities it is the age
#' at the time of capture/recapture.
#' 
#' \code{groups} is a vector of variable names that are contained in
#' \code{data}.  Each must be a factor variable. A group is created for each
#' unique combination of the levels of the factor variables.  In the first
#' example given below \code{groups=c("sex","age","region")}. which creates
#' groups defined by the levels of \code{sex}, \code{age} and \code{region}.
#' There should be 2(sexes)*3(ages)*4(regions)=24 groups but in actuality there
#' are only 16 in the data because there are only 2 age groups for each sex.
#' Age group 1 and 2 for M and age groups 2 and 3 for F.  This was done to
#' demonstrate that the code will only use groups that have 1 or more capture
#' histories unless \code{allgroups=TRUE}.
#' 
#' The argument \code{age.var=2} specifies that the second grouping variable in
#' \code{groups} represents an age variable.  It could have been named
#' something different than age. If a variable in \code{groups} is named
#' \code{age} then it is not necessary to specify \code{age.var}.
#' \code{initial.age} specifies that the age at first capture of the age levels
#' is 0,1 and 2 while the age classes were designated as 1,2,3. The actual ages
#' for the age classes do not have to be sequential or ordered, but ordering
#' will cause less confusion.  Thus levels 1,2,3 could represent initial ages
#' of 0,4,6 or 6,0,4. The argument age.unit is the amount an animal ages for
#' each unit of time and the default is 1.  The default for \code{initial.age}
#' is 0 for each group, in which case, \code{age} represents time since marking
#' (first capture) rather than the actual age of the animal.
#' 
#' @param data A data frame with at least one field named \code{ch} which is
#' the capture (encounter) history stored as a character string. \code{data}
#' can also have a field \code{freq} which is the number of animals with that
#' capture history. The default structure is freq=1 and it need not be included
#' in the dataframe. \code{data} can also contain an arbitrary number of
#' covariates specific to animals with that capture history.
#' @param begin.time Time of first capture occasion or vector of times if
#' different for each group
#' @param model Type of analysis model. See \code{\link{mark}} for a list of
#' possible values for \code{model}
#' @param mixtures Number of mixtures in closed capture models with
#' heterogeneity or number of secondary samples for MultScaleOcc model
#' @param groups Vector of factor variable names (in double quotes) in
#' \code{data} that will be used to create groups in the data. A group is
#' created for each unique combination of the levels of the factor variables in
#' the list.
#' @param allgroups Logical variable; if TRUE, all groups are created from
#' factors defined in \code{groups} even if there are no observations in the
#' group
#' @param age.var An index in vector \code{groups} for a variable (if any) for
#' age
#' @param initial.ages A vector of initial ages that contains a value for each
#' level of the age variable \code{groups[age.var]}
#' @param age.unit Increment of age for each increment of time as defined by
#' \code{time.intervals}
#' @param time.intervals Vector of lengths of time between capture occasions
#' @param nocc number of occasions for Nest type; either nocc or time.intervals
#' must be specified
#' @param strata.labels vector of single character values used in capture
#' history(ch) for ORDMS, CRDMS, RDMSOccRepro models; it can contain one more value beyond what is
#' in ch for an unobservable state except for RDMSOccRepro which is used to specify strata ordering (eg 0 not-occupied, 1 occupied no repro, 2 occupied with repro.
#' @param counts named list of numeric vectors (one group) or matrices (>1
#' group) containing counts for mark-resight models
#' @param reverse if set to TRUE, will reverse timing of transition (Psi) and
#' survival (S) in Multistratum models
#' @return processed.data (a list with the following elements)
#' \item{data}{original raw dataframe with group factor variable added if
#' groups were defined} \item{model}{type of analysis model (eg, "CJS",
#' "Burnham", "Barker")} \item{freq}{a dataframe of frequencies (same number of
#' rows as data, number of columns is the number of groups in the data. The
#' column names are the group labels representing the unique groups that have
#' one or more capture histories.} \item{nocc}{number of capture occasions}
#' \item{time.intervals}{length of time intervals between capture occasions}
#' \item{begin.time}{time of first capture occasion} \item{age.unit}{increment
#' of age for each increment of time} \item{initial.ages}{an initial age for
#' each group in the data; Note that this is not the original argument but is a
#' vector with the initial age for each group. In the first example below
#' \code{proc.example.data$initial.ages} is a vector with 16 elements as
#' follows 0 1 1 2 0 1 1 2 0 1 1 2 0 1 1 2} \item{nstrata}{number of strata in
#' Multistrata models} \item{strata.labels}{vector of alphabetic characters
#' used to identify strata in Multistrata models}
#' \item{group.covariates}{factor covariates used to define groups}
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{import.chdata}}, \code{\link{dipper}},
#' \code{\link{edwards.eberhardt}}, \code{\link{example.data}}
#' @keywords utility
#' @examples
#' 
#' data(example.data)
#' proc.example.data=process.data(data=example.data,begin.time=1980,
#' groups=c("sex","age","region"),
#' age.var=2,initial.age=c(0,1,2))
#' 
#' data(dipper)
#' dipper.process=process.data(dipper)
#' 
process.data <-
function(data,begin.time=1,model="CJS",mixtures=1,groups=NULL,allgroups=FALSE,age.var=NULL,
initial.ages=c(0),age.unit=1,time.intervals=NULL,nocc=NULL,strata.labels=NULL,counts=NULL,reverse=FALSE)
{
robust.occasions<-function(times)
{
   if(times[1] !=0 | times[length(times)]!=0)
      stop("\nIncorrect structure for time intervals with robust design. Time intervals must begin and end with a zero.\n")
   merge_int=unlist(strsplit(paste(as.numeric(times>0),collapse=""),"0"))
   if(length(merge_int[merge_int!=""])!=length(times[times>0]))
	   stop("\nIncorrect structure for time intervals with robust design. Must have at least 2 secondary occasions per primary period.\n")
   times.vec=match(as.character(times),"0")
   nocc.secondary=nchar(unlist(strsplit(paste(as.character(times.vec),collapse=""),"NA")))
   nocc=length(nocc.secondary)
   nocc.secondary=nocc.secondary+1
   if(any(nocc.secondary==0))
      stop("\nIncorrect setup for time intervals. Must have at least one zero between non-zero intervals.\n")
   return(list(nocc=nocc,nocc.secondary=nocc.secondary))
}
   dataname=substitute(data)
#
#  handle Multistrata reversal if needed
#
   if(reverse)
   {
	   if(model!="Multistrata")
		   stop("reverse can only be set TRUE for the Multistrata model")
	   else
	   {
		   if(is.null(time.intervals))time.intervals=rep(1,nchar(data$ch[1])-1)
		   data$ch=sapply(strsplit(data$ch,""),paste,collapse="0") 
		   time.intervals=as.vector(matrix(c(rep(0,length(time.intervals)),time.intervals),ncol=length(time.intervals),nrow=2,byrow=TRUE))
	   }
   }
#
#   If model="Nest" a completely different structure is used.  No ch is used; fate used instead
#
   if(model=="Nest")
   {
      begin.time=1
      nstrata=1
      strata.labels=""
      if(!all(c("FirstFound","LastPresent","LastChecked","Fate") %in% names(data)))
         stop("data should contain fields: FirstFound, LastPresent, LastChecked, Fate. One or more are missing.")
      if(is.null(time.intervals))
      {
         if(is.null(nocc)) stop("nocc or time.intervals must be set for Nest model")
         time.intervals=rep(1,nocc)
      }
      else
      {
         if(is.null(nocc)) nocc=length(time.intervals)
         if(length(time.intervals)!=nocc)stop("length of time intervals must match nocc for Nest model")
      }
      if(any(data$FirstFound>nocc | data$FirstFound< 1))
         stop("One or more FirstFound values greater than number of occasions or less than 1")
      if(any(data$LastChecked>nocc| data$LastChecked< 1))
         stop("One or more LastChecked values greater than number of occasions or less than 1")
      if(any(data$LastPresent>nocc | data$LastPresent<1))
         stop("One or more LastPresent values greater than number of occasions or less than 1")
   }
   else
   {
#
#  Compute number of occasions and check validity of model
#
      if(is.null(data$ch))
       stop("Field ch is missing in ",substitute(data))
      ch.lengths=nchar(data$ch)
      nocc=median(ch.lengths)
      if(any(ch.lengths!=nocc))
      {
           stop(paste("\nCapture history length is not constant. ch must be a character string",
               "\n row numbers with incorrect ch length",paste(row.names(data[ch.lengths!=nocc,]),collapse=","),"\n"))
      }
  }
#
#  Setup model
#
   model.list=setup.model(model,nocc,mixtures)
#
#  data checks
#
   if(model!="Nest")
   {
	   if(!is.character(data$ch))
	   {
		   if(is.factor(data$ch))
			   stop("\nch field is a factor and must be a character string\n")
		   else
			   stop("\nch field must be a character string\n")
	   } else
	   if(!model.list$occupancy & !model.list$model %in% c("LogitNormalMR","IELogitNormalMR"))
		   if(any(sapply(strsplit(data$ch,""),function(x) all(x=="0"))))
			   stop("\nall 0 ch encountered. MARK will not accept them\n")
	   if(!is.null(data$freq))
	   {
		   if(!is.numeric(data$freq))
			   stop("\n freq field must be numeric\n")
	   }  
   }
   #
#  If multistrata design, determine number of strata and their labels
#  Make sure multistrata designs have at least 2 strata
#
   if(model!="Nest")
   {
      ch.values=unique(unlist(strsplit(data$ch,"")))
      if(model=="MSLiveDead")
         inp.strata.labels=sort(ch.values[!(ch.values %in% c("0",".","1"))])
      else
         inp.strata.labels=sort(ch.values[!(ch.values %in% c("0","."))])
	  if(model%in%c("RDMSOpenMisClass","RDMSMisClass","RDMS2MisClass","RDMSOpenMCSeas"))
		  inp.strata.labels=inp.strata.labels[!inp.strata.labels%in%"u"]
      nstrata = length(inp.strata.labels)                  
      if(model.list$strata)
      {
        if(is.null(strata.labels)) 
        {
           strata.labels=inp.strata.labels
        }
        else
        {
		   if(any(!nchar(strata.labels)==1))stop("strata.labels must be a single character")
           nstrata=length(strata.labels)
           if(!all(inp.strata.labels %in% strata.labels))
              stop(paste("Some strata labels in data",paste(inp.strata.labels,collapse=","),"are not in strata.labels"))
           if(sum(as.numeric(strata.labels %in% inp.strata.labels))< (nstrata-1))
              message("Note: More than one non-observable state has been specified")
        }
        if(nstrata<2)stop("\nAny multistrata model must have at least 2 strata\n")
      } else
      {              
         nstrata=1
         if(!is.null(model.list$occupancy) && model.list$occupancy)
         {
            if(model == "OccupRPoisson" | model=="OccupRNegBin")
            {
               if(any(!ch.values%in%c(0:9,".")))
                  stop(paste("\nIncorrect count values in data:",paste(ch.values,collapse=""),"\n",sep=""))
            }
            else
            {
               if(any(!ch.values%in%c(".","0","1","2")))
                    stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
            }
         }
         else
         {
            if(any(!ch.values%in%c("0","1",".")))
            {
               if(substr(model,1,6)!="RDBark" & model!="RDBarkHug" & model!="Barker" & model!="MSOccupancy" &model!="PoissonMR"&model!="UnIdPoissonMR")
                  stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
               else
				  if(model%in%c("PoissonMR","UnIdPoissonMR"))
				  {			  
					  if(any(!ch.values%in%c(".",as.character(0:9),"+","-")))
						  stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
				  } else
				  {
				      if(any(!ch.values%in%c(".","0","1","2")))
                         stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
			      }
            }
         }
      }
   }
#
#  If this is a robust design, compute number of primary and secondary occasions
#
   if(model.list$robust &model!="MultScalOcc")
   {
      if(is.null(time.intervals))
         stop("\nTime intervals must be specified for a robust design\n")
      else
      {
		 if(substr(model,1,6)=="RDBark")
			  nocc.list=robust.occasions(time.intervals[-length(time.intervals)])
		 else
		   nocc.list=robust.occasions(time.intervals)
		 nocc=nocc.list$nocc
         nocc.secondary=nocc.list$nocc.secondary
         if(any(nchar(data$ch)/model.list$divisor !=sum(nocc.secondary)))
             stop("Incorrect number of time intervals. One or more capture history lengths do not match time interval structure.")

      }
      num=model.list$num
   }
   else
   {
	   nocc=model.list$nocc
	   nocc.secondary=NULL
	   num=model.list$num
	   if(model=="MultScalOcc"){
		   nocc.secondary=mixtures
		   if(nocc==nocc.secondary*(nocc %/% nocc.secondary))
		   {
		      nocc=nocc/nocc.secondary
			  nocc.secondary=rep(nocc.secondary,nocc)
		   }
		   else
			  stop("# of mixtures (secondary samples) not an even multiple of ch length")
		   model.list$mixtures=1
	   }
#
#     If time intervals specified make sure there are nocc-1 of them
#     If none specified assume they are 1
#
      if(is.null(time.intervals))
          if(model=="MultScalOcc")
			  time.intervals=rep(1,nocc.secondary[1]*nocc-1)
	      else
		      time.intervals=rep(1,nocc+model.list$num)
      else
         if(length(time.intervals)!=(nocc+num))
             stop("Incorrect number of time intervals")
   }
   mixtures=model.list$mixtures
#
#  Get number of factors to create groups
#
   if(is.null(groups))
     number.of.factors=0
   else
     number.of.factors=length(groups)
#
#  Get number of records in data set
#
number.of.ch=dim(data)[1]
#
#  See if there is already a freq variable in data set
#
if(!is.null(data$Freq)) names(data)[which("Freq"== names(data))]="freq"
has.freq=!is.null(data$freq)
#
#  If there are no factors then
#     if already has freq variable return the input data set as a list
#     otherwise add the freq variable with each value = 1 and return as a list
#  If model=js, then add dummy data for non-captured 
#
if(number.of.factors==0)
{
   if(has.freq)
   {
#       if(model=="js")
#       {
#          data=add.dummy.data(data,nocc=nocc,group.covariates=NULL)     
#          number.of.ch=dim(data)[1]
#       }
       return(list(data=data,model=model,mixtures=mixtures,
                   freq=matrix(data$freq,ncol=1,dimnames=list(1:number.of.ch,"group1")),
                   nocc=nocc, nocc.secondary=nocc.secondary,time.intervals=time.intervals,begin.time=begin.time,
                   age.unit=age.unit,initial.ages=initial.ages[1],group.covariates=NULL,nstrata=nstrata,
                   strata.labels=strata.labels,counts=counts,reverse=reverse))
   }
   else
   {
       data$freq=rep(1,number.of.ch)
#       if(model=="js")
#       {
#          data=add.dummy.data(data,nocc=nocc,group.covariates=NULL)            
#          number.of.ch=dim(data)[1]
#       }
       return(list(data=data,model=model,mixtures=mixtures,
                   freq=matrix(rep(1,number.of.ch),ncol=1,dimnames=list(1:number.of.ch,"group1")),
                   nocc=nocc,  nocc.secondary=nocc.secondary, time.intervals=time.intervals,begin.time=begin.time,
                   age.unit=age.unit,initial.ages=initial.ages[1],group.covariates=NULL,nstrata=nstrata,
                   strata.labels=strata.labels,counts=counts,reverse=reverse))
   }
}
#
#   If there are one or more in the group factor list then
#     make sure each is a factor variable in the data set and compute number
#         of levels for each factor, cumlevels and factor matrix
#     if not a factor variable - stop with error message
# 
else
{
  if(!has.freq)
       data$freq=rep(1,number.of.ch)
  number.of.groups=1
#  var=rep(" ",number.of.factors)
  n.levels=rep(0,number.of.factors)
  facmat=NULL
  faclabs=list()
  for (i in 1:number.of.factors)
  {
    vari=data[,groups[i]]
#    var[i]=paste(dataname,"$",groups[i],sep="")
#    vari=eval.parent(parse(text=var[i]))
    if(!is.factor(vari))
	{
		warning(paste("\n ",groups[i]," is not a factor variable. Coercing to factor.\n"))
	    data[,groups[i]]=factor(data[,groups[i]])
		vari=data[,groups[i]]
	}
    n.levels[i]=length(levels(vari))
    facmat=cbind(facmat,as.numeric(vari)-1)
    faclabs[[i]]=levels(vari)
  }
  cumlevels=cumprod(n.levels)
  number.of.groups=cumlevels[length(cumlevels)]

#  If age.var is specified, make sure it is valid and that the number of 
#  initial.ages matches number of levels of identified variable
#
   if(is.null(age.var))
      age.var=match("age",groups)
   if(!is.na(age.var))
   {
      if(age.var>length(groups) | age.var<1)
         stop("Invalid age variable. Must be between 1 and ",length(groups))
      if(is.null(initial.ages))
         stop("initial.ages must be specified if age.var is specified")
      else
         if(!is.numeric(initial.ages) | (length(initial.ages)!=n.levels[age.var] &length(initial.ages)>1) )
           stop(paste("intial.ages must be numeric and match length of levels of",groups[age.var]))
      if(age.unit<=0 | !is.numeric(age.unit)) stop("age.unit must be numeric and >0")
   }
   else
	   if(!is.null(initial.ages)&&length(initial.ages)>1)
		   stop("initial.ages vector specified but no age variable was identified in groups; see age.var")
#
#  Next compute the group number for each capture history
#
   if(number.of.factors==1)
      data$group=facmat+1
   else
      if(number.of.factors==2)
         data$group=facmat[,2]*cumlevels[1]+facmat[,1]+1
      else
         data$group=facmat[,2:number.of.factors]%*%cumlevels[1:(number.of.factors-1)]+facmat[,1]+1
#
#  Next create frequency matrix for groups   
#
  freqmat=matrix(0,nrow=number.of.ch,ncol=number.of.groups)
  for(i in 1:number.of.ch)
  {
     freqmat[i,data$group[i]]=data$freq[i]
  }
#
#  If allgroups=FALSE, recompute number of groups and group number based on groups with 1 or more capture histories
#
  if(!allgroups)
  {
     test.freq=freqmat
     test.freq[test.freq!=0]=1
     chcounts = apply(test.freq, 2, sum)
     newgroups=rep(0,number.of.groups)
     index=1
     for (i in 1:number.of.groups)
        if(chcounts[i]>0)
        {
           newgroups[i]=index
           index=index+1
        }     
     data$group=as.factor(newgroups[data$group])
     freqmat=freqmat[,chcounts>0]
     number.of.groups=index-1
  }
#
#  Check to make sure length of begin.time is either 1 or equal to the
#  number of groups
#
  if(length(begin.time)!=1 & length(begin.time)!=number.of.groups)
    stop("length of begin.time must either be 1 or match number of groups")
#
#  Create group labels
#  
  labs=expand.grid(faclabs)
  if(!allgroups)labs=as.matrix(labs[chcounts>0,])
#
#  If age.var has not been set, initial ages are set to 0
#
  if(is.na(age.var))
    init.ages=rep(initial.ages[1],number.of.groups)
  else
  {
    if(length(initial.ages)==1)
       initial.ages=rep(initial.ages,length(levels(as.factor(labs[,age.var]))))
    init.ages = initial.ages[as.numeric(factor(labs[,age.var],levels=unique(faclabs[[age.var]])))]
  }
  grouplabs=rep(" ",number.of.groups)
  for (i in 1:number.of.groups)
     grouplabs[i]=paste(groups,labs[i,],sep="",collapse=".") 
  freqmat=as.data.frame(freqmat)
  names(freqmat)=grouplabs
#
#  Store labs as group covariates; set levels to the same as in data
#  
  group.covariates=as.data.frame(labs)
  names(group.covariates)=groups
  for (i in 1:dim(group.covariates)[2])
     group.covariates[,i]=factor(group.covariates[,i],levels=levels(data[,groups[i]]))
#
# Return data as a list with original dataframe and frequency matrix
#
#  if(model=="js")
#     data=add.dummy.data(data,nocc,group.covariates)     
#  else
  if(!has.freq)data$freq=NULL
  return(list(data=data,model=model,mixtures=mixtures,freq=freqmat,
                   nocc=nocc, nocc.secondary=nocc.secondary, time.intervals=time.intervals,begin.time=begin.time,
                   age.unit=age.unit,initial.ages=init.ages,
                   group.covariates=group.covariates,nstrata=nstrata,
                   strata.labels=strata.labels,counts=counts,reverse=reverse))
}
}
