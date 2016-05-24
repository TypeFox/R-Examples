#' Process encounter history dataframe for MARK analysis
#' 
#' Prior to analyzing the data, this function initializes several variables
#' (e.g., number of capture occasions, time intervals) that are often specific
#' to the capture-recapture model being fitted to the data.  It also is used to
#' 1) define groups in the data that represent different levels of one or morestrata.labels
#' factor covariates (e.g., sex), 2) define time intervals between capture
#' occasions (if not 1), and 3) create an age structure for the data, if any.
#' 
#' For examples of \code{data}, see \code{\link{dipper}}. The structure of the
#' encounter history and the analysis depends on the analysis model to some
#' extent. Thus, it is necessary to process a dataframe with the encounter
#' history (\code{ch}) and a chosen \code{model} to define the relevant values.
#' For example, number of capture occasions (\code{nocc}) is automatically
#' computed based on the length of the encounter history (\code{ch}) in
#' \code{data}. Currently, only 2 types of models are accepted in marked: cjs and js.  
#' The default time interval is unit time (1) and if this is
#' adequate, the function will assign the appropriate length.  A processed data
#' frame can only be analyzed using the model that was specified.  The
#' \code{model} value is used by the functions \code{\link{make.design.data}}
#' and \code{\link{crm}} to define the model structure as it relates to the
#' data. Thus, if the data are going to be analysed with different underlying
#' models, create different processed data sets with the model name as an
#' extension.  For example, \code{dipper.cjs=process.data(dipper)}.
#' 
#' This function will report inconsistencies in the lengths of the capture
#' history values and when invalid entries are given in the capture history.
#' 
#' The argument \code{begin.time} specifies the time for the first capture
#' occasion and not the first time the particular animal was caught or releaed.  This is used in creating the levels of the time factor variable
#' in the design data and for labelling parameters. If the \code{begin.time}
#' varies by group, enter a vector of times with one for each group. It will add a field
#' begin.time to the data with the value for each individual.  You can also specify a
#' begin.time field in the data allowing each animal to have a unique begin.time. Note that
#' the time values for survivals are based on the beginning of the survival
#' interval and capture probabilities are labeled based on the time of the
#' capture occasion.  Likewise, age labels for survival are the ages at the
#' beginning times of the intervals and for capture probabilities it is the age
#' at the time of capture/recapture. 
#' 
#' The time.intervals argument can either be a vector of lengths of times for each interval between occasions
#' that is constant for all animals or a matrix which has a row for each animal and a column for each
#' interval which lets the intervals vary by animals. These intervals are used to construct the design data
#' and are used for the field time.interval which is used to adjust parameters like Phi and S to a constant per
#' unit time interval (eg annual survival rates). On occasion it can be useful to leave the time.interval to 
#' remain at default of 1 or some other vector of time.intervals to construct the design data and then modify the 
#' time.interval value in the design data.  For example, assume that cohort marking and release is done between
#' sampling occasions.  The initial survival from release to the next sampling occasion may vary by release
#' cohort, but the remainder of the surivivals are between sampling occasions.  In that case it is easier to 
#' let time.interval=1 (assuming unit interval (eg year) between sampling occasions but then modifying ddl$Phi$time.interval
#' to the value for the first interval after each release to be the partial year from release to next sampling occasion. In
#' this way everything is labelled with annual quantities but the first partial year survival is adjusted to an annual rate.
#' 
#' Note that if you specify time.intervals as a matrix, then accumulate is set to FALSE so that the number of 
#' rows in the data can be checked against the number of rows in the time.intervals matrix and thus data cannot be
#' accumulated because at present it doesn't use values of time.intervals to determine which records can be accumulated.
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
#' \code{initial.age} specifies that the age at first capture of the age levels. For example
#' initial.age=0:2 specifies that the initial.ages are 0,1 and 2 for the age class levels 
#'  designated as 1,2,3. The actual ages
#' for the age classes do not have to be sequential or ordered, but ordering
#' will cause less confusion.  Thus levels 1,2,3 could represent initial ages
#' of 0,4,6 or 6,0,4. The default for \code{initial.age}
#' is 0 for each group, in which case, \code{age} represents time since marking
#' (first capture) rather than the actual age of the animal. If the data contains an initial.age field
#' then it overrides any other values and lets each animal have a unique initial.age at first capture/release.
#' 
#' The following variable names are reserved and should be used as follows:
#' id (animal id)
#' ch(capture history)
#' freq (number of animals with that ch/data)
#' The following variable names are reserved and should not be used in the data:
#' occ,age,time,cohort,Age,Time,Cohort,Y,Z,initial.age,begin.time,time.interval,fix
#' 
#' @aliases process.data accumulate_data
#' @usage 	process.data(data,begin.time=1,model="CJS",mixtures=1,groups=NULL,
#'                        allgroups=FALSE,age.var=NULL,initial.ages=c(0), 
#'                        time.intervals=NULL,nocc=NULL,accumulate=TRUE,
#'                        strata.labels=NULL)
#' 
#'          accumulate_data(data)
#' 
#' @param data A data frame with at least one field named \code{ch} which is
#' the capture (encounter) history stored as a character string. \code{data}
#' can also have a field \code{freq} which is the number of animals with that
#' capture history. The default structure is freq=1 and it need not be included
#' in the dataframe. \code{data} can also contain an arbitrary number of
#' covariates specific to animals with that capture history.
#' @param begin.time Time of first capture occasion or vector of times if
#' different for each group
#' @param model Type of analysis model. 
#' @param mixtures Number of mixtures in closed capture models with
#' heterogeneity
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
#' level of the age variable \code{groups[age.var]}; if the data contain an initial.age field then it will be used instead.
#' @param time.intervals Vector of lengths of time between capture occasions or matrix of time intervals with a row for each animal and column for each interval between occasions.
#' @param nocc number of occasions for Nest type; either nocc or time.intervals
#' must be specified
#' @param accumulate if TRUE, aggregates data with same values and creates freq field for count of records
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata; for any unobserved strata p=0.
#' @return from \code{process.data} processed.data (a list with the following elements)
#' \item{data}{original raw dataframe with group factor variable added if
#' groups were defined} \item{model}{type of analysis model (eg, "cjs" or "js")}
#' \item{freq}{a dataframe of frequencies (same # of rows
#' as data, number of columns is the number of groups in the data. The column
#' names are the group labels representing the unique groups that have one or
#' more capture histories.} \item{nocc}{number of capture occasions}
#' \item{time.intervals}{length of time intervals between capture occasions}
#' \item{begin.time}{time of first capture occasion} \item{initial.ages}{an initial age for
#' each group in the data; Note that this is not the original argument but is a
#' vector with the initial age for each group. In the first example below
#' \code{proc.example.data$initial.ages} is a vector with 16 elements as
#' follows 0 1 1 2 0 1 1 2 0 1 1 2 0 1 1 2} \item{group.covariates}{factor covariates used to define groups}
#' from accumulate_data a dataframe with same column structure as argument with addition of freq (if not any)
#' and reduced to unique rows with freq accumulating number of records. 
#' @author Jeff Laake
#' @export process.data accumulate_data
#' @seealso \code{\link{dipper}},\code{\link{crm}}
#' @keywords utility
#' @examples
#' 
#' 
#' data(dipper)
#' dipper.process=process.data(dipper)
#' accumulate_data(dipper)
#' 
process.data <-
function(data,begin.time=1,model="CJS",mixtures=1,groups=NULL,allgroups=FALSE,age.var=NULL,
initial.ages=c(0),time.intervals=NULL,nocc=NULL,accumulate=TRUE,strata.labels=NULL)
{
   model=toupper(model)
   dataname=substitute(data)
  #
  #  Compute number of occasions and check validity of model
  #
   if(is.null(data$ch))
     stop("Field ch is missing in ",substitute(data))
   # If ch is not comma delimited, turn into comma delimited
   if(length(grep(",",data$ch[1]))==0)
	   data$ch=sapply(strsplit(data$ch,""),paste,collapse=",")
#   data$ch=toupper(data$ch)
   ch.lengths=sapply(strsplit(data$ch,","),length)
   nocc=median(ch.lengths)
   if(any(ch.lengths!=nocc))
   {
        stop(paste("\nCapture history length is not constant. ch must be a character string",
            "\n row numbers with incorrect ch length",paste(row.names(data[ch.lengths!=nocc,]),collapse=","),"\n"))
   }
   #
   #  Setup model
   #
   model.list=setup.model(model,nocc,mixtures)
   # regardless of user input for accumulate if it is FALSE in model.list set to FALSE;
   # the bayesian models cannot deal with accumulation
   if(!model.list$accumulate)accumulate=FALSE
   ch.values=unique(unlist(strsplit(data$ch,",")))
   nocc=model.list$nocc
   nocc.secondary=NULL
   num=model.list$num
   if(model.list$strata&is.null(strata.labels)&!model%in%c("HMMCJS2TL","HMMCJS1TL"))
	   stop("\nstrata.labels must be specified for stratified models\n")
   if(model.list$IShmm)
   {
	   model.list=setupHMM(model.list,model,strata.labels)
	   if(model.list$strata)strata.labels=model.list$hmm$strata.labels
   }
   #  If no strata in model then only 0,1 are acceptable values
   if(!model.list$strata)
   {
	   if(any(!ch.values%in%c("0","1")))
		   stop(paste("\nIncorrect ch values in data:",paste(ch.values,collapse=""),"\n",sep=""))
   } else
   {
	   # Get unique ch values and use as strata.labels unless they are specified   
	   inp.strata.labels=sort(ch.values[!(ch.values %in% c("0"))])
	   if(substr(model,1,4)%in%c("HMMU","MVMS"))
	   {
		   if(substr(model,1,4)=="MVMS")
			   uindex=grep("u",inp.strata.labels)
		   else
			   uindex=grep("U",inp.strata.labels)
		   if(length(uindex)>0)inp.strata.labels=inp.strata.labels[-uindex]
	   }
	   nstrata = length(inp.strata.labels)                  
	   if(is.null(strata.labels))
		   strata.labels=inp.strata.labels
	   # If strata specified as a vector test otherwise pass it through
	   if(is.vector(strata.labels)&!model%in%c("HMMCJS2TL","HMMCJS1TL"))
	   {
		   strata.labels=c(strata.labels[strata.labels%in%inp.strata.labels],strata.labels[!strata.labels%in%inp.strata.labels])
		   nstrata=length(strata.labels)
		   unobserved=nstrata-length(inp.strata.labels)
		   if(unobserved<0 | !all(inp.strata.labels%in%strata.labels))
			   stop(paste("\nSome of the strata in the data\n",paste(inp.strata.labels,collapse=","),"\nnot specified in strata.labels\n",paste(strata.labels,collapse=","),"\n"))
		   if(nstrata<2)stop("\nAny multistrata model must have at least 2 strata\n")
	   } else
		   unobserved=0
   }
   #
   #     If time intervals specified make sure there are nocc-1 of them if a vector
   #     and if a matrix rows must match number of animals and # cols = nocc+num
   #     If none specified assume they are 1
   #
   if(is.null(time.intervals))
      time.intervals=rep(1,nocc+model.list$num)
   else
	   if(is.vector(time.intervals))
	   {
		   if(length(time.intervals)!=(nocc+num))
			   stop("Incorrect number of time intervals")
	   }else
	   {
		   if(is.matrix(time.intervals))
		   {
			   if(ncol(time.intervals)!=(nocc+num))
				   stop(paste("Incorrect number of columns in time.intervals. Should be:",nocc+num))
			   if(nrow(time.intervals)!=nrow(data))
				   stop(paste("Incorrect number of rows in time.intervals. Should be:",nrow(data)))
		   }else
			   stop("\ntime.intervals must be either a vector or matrix")
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
   # Accumulate non-unique data; if accumulate=F and dataframe has freq>1 expand and add id field
   # if null; if cjs type remove all histories that are initiated on last occasion
   # 
   if(model.list$cjs)
   {
	  chproc=process.ch(data$ch,freq=1,all=FALSE)
	  if(!is.matrix(time.intervals))data=data[chproc$first!=chproc$nocc,,drop=FALSE]
   }	
   if(!is.null(data$Freq)) names(data)[which("Freq"== names(data))]="freq"
   if(accumulate&is.matrix(time.intervals))
   {
	   accumulate=FALSE   
	   message("ignoring accumulate=TRUE value because time.intervals specified as a matrix")
   }
   if(accumulate)
    	data=accumulate_data(data)
   else
   {
	   if(is.null(data$freq))
	    	data$freq=1
	   else
	   {
	    	data=data[rep(1:nrow(data),times=data$freq),]
		    data$freq=1
	   }
   }
   if(is.null(data$id))
    	data$id=factor(1:nrow(data))
   else
        stop("data argument cannot contain a variable named id. It is reserved.")	
   #
   #  Get number of records in data set
   #
   number.of.ch=nrow(data)
   # start - for each ch, the first non-zero x value and the occasion of the first non-zero value
   #  if strata.labels then the first non-zreo x is matched to strata.labels to exclude possible "U" values
   if(is.null(model.list$hmm$strata.labels))
       start=t(sapply(data$ch,function(x){
					xx=strsplit(x,",")[[1]]
					ich=min(which(strsplit(x,",")[[1]]!="0"))
					return(c(as.numeric(factor(xx[ich],levels=model.list$hmm$ObsLevels))-1,ich))
				}))
   else
	   if(is.null(model.list$hmm$obs_strata_map))
	       start=t(sapply(data$ch,function(x){
						   xx=strsplit(x,",")[[1]]
						   ich=min(which(strsplit(x,",")[[1]]!="0"))
						   return(c(match(xx[ich],model.list$hmm$strata.labels),ich))
					   }))
       else
		   start=t(sapply(data$ch,function(x){
							   xx=strsplit(x,",")[[1]]
							   ich=min(which(strsplit(x,",")[[1]]!="0"))
							   return(c(model.list$hmm$obs_strata_map[match(xx[ich],model.list$hmm$ObsLevels)],ich))
						   }))
   
   # create encounter history matrix
   if(model.list$IShmm)
      ehmat=t(sapply(strsplit(data$ch,","),function(x) as.numeric(factor(x,levels=model.list$hmm$ObsLevels))))
   else
	  ehmat=process.ch(data$ch)$chmat
   #
   #  If there are no factors then
   #     if already has freq variable return the input data set as a list
   #     otherwise add the freq variable with each value = 1 and return as a list
   #  If model=js, then add dummy data for non-captured 
   #
   if(number.of.factors==0)
   {
       if(model=="JS")
       {
          data=add.dummy.data(data,nocc=nocc,group.covariates=NULL)     
          number.of.ch=nrow(data)
		  data$id=factor(1:nrow(data))
       }
	   if(length(begin.time)>1)stop("\nbegin.time has more than one value and no groups were specified")
       plist=list(data=data,model=model,mixtures=mixtures,
                   freq=matrix(data$freq,ncol=1,dimnames=list(1:number.of.ch,"group1")),
                   nocc=nocc, nocc.secondary=nocc.secondary,time.intervals=time.intervals,begin.time=begin.time,
                   initial.ages=initial.ages[1],group.covariates=NULL,start=start,ehmat=ehmat)
 	   if(model.list$IShmm)
		   plist=c(plist,model.list$hmm)
	   else
	       if(model.list$strata)
			  plist=c(plist,list(strata=model.list$strata,strata.labels=strata.labels,unobserved=unobserved))
	   return(plist)
   }
   #
   #   If there are one or more in the group factor list then
   #     make sure each is a factor variable in the data set and compute number
   #         of levels for each factor, cumlevels and factor matrix
   #     if not a factor variable - stop with error message
   # 
   else
   {
      number.of.groups=1
      n.levels=rep(0,number.of.factors)
      facmat=NULL
      faclabs=list()
      for (i in 1:number.of.factors)
      {
        vari=data[,groups[i]]
        if(!is.factor(vari))
            stop(paste("\n ",groups[i]," is not a factor variable\n"))
        else
        {
            n.levels[i]=length(levels(vari))
            facmat=cbind(facmat,as.numeric(vari)-1)
            faclabs[[i]]=levels(vari)
        }       
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
      }
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
         freqmat[i,data$group[i]]=data$freq[i]
      #
      #  If allgroups=FALSE, recompute number of groups and group number based on groups with 1 or more capture histories
      #
      if(!allgroups)
      {
        test.freq=freqmat
        test.freq[test.freq!=0]=1
        counts = apply(test.freq, 2, sum)
        newgroups=rep(0,number.of.groups)
        index=1
        for (i in 1:number.of.groups)
           if(counts[i]>0)
           {
              newgroups[i]=index
              index=index+1
           }     
        data$group=as.factor(newgroups[data$group])
        freqmat=freqmat[,counts>0]
        number.of.groups=index-1
      } 
      #
      #  Create group labels
      #  
      labs=expand.grid(faclabs)
      if(!allgroups)labs=as.matrix(labs[counts>0,])
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
        if(model=="JS")
          data=add.dummy.data(data,nocc,group.covariates)     
	    data$id=factor(1:nrow(data))
	    if(is.null(data$initial.age)) data$initial.age=init.ages[data$group]
		#
		#  Check to make sure length of begin.time is either 1 or equal to the
		#  number of groups
		#
		if(length(begin.time)!=1)
		if(length(begin.time)!=number.of.groups)
			stop("length of begin.time must either be 1 or match number of groups")
		else
		    data$begin.time=begin.time[data$group]
		plist=list(data=data,model=model,mixtures=mixtures,freq=freqmat,
				nocc=nocc, nocc.secondary=nocc.secondary, time.intervals=time.intervals,begin.time=begin.time,
				initial.ages=init.ages,group.covariates=group.covariates,start=start,ehmat=ehmat)
        if(model.list$strata)plist=c(plist,list(strata=model.list$strata,strata.labels=model.list$strata.labels,unobserved=unobserved))
        if(!is.null(model.list$hmm)) 
		{
			if(!is.null(model.list$hmm$strata.labels))
		    {
		        plist$strata.labels=model.list$hmm$strata.labels	
			    plist=c(plist,model.list$hmm[!names(model.list$hmm)%in%"strata.labels"])
			}
			else
				plist=c(plist,model.list$hmm)
		}
        return(plist) 
    }
}
add.dummy.data=function(data,nocc,group.covariates)
{
	if(is.null(group.covariates))
		number.of.groups=1
	else
		number.of.groups=nrow(group.covariates)
	if(!is.null(group.covariates))
	{
		xlist=split(data,data[,names(group.covariates),drop=FALSE])
		xlist=xlist[as.vector(sapply(xlist,function(x) nrow(x)))>0]
	} else
	{                              
		xlist=data
	}
	numvar=sapply(data,is.numeric)
	numvar=numvar[names(data)!="freq"]
	if(any(numvar))
	{
		numvar=names(data[,names(data)!="freq"])[numvar]
		if(!is.null(group.covariates))
		{
			xmeans=sapply(xlist,function(x) sapply(subset(x,select=numvar),mean))
			if(length(numvar)==1)
			{
				dd=data.frame(group=1:nrow(group.covariates))
				dd[,numvar]=xmeans
			}else
			{
				dd=data.frame(group=1:nrow(group.covariates))
				dd=cbind(dd,t(xmeans))
			}
			xx=merge(cbind(data.frame(group=1:nrow(group.covariates),group.covariates)),
					dd,by.x="group",all.x=TRUE)
		}
		else
		{
			xmeans=colMeans(subset(data,select=numvar))    
			if(ncol(t(xmeans))==0)
				xx=data.frame(N=1)
			else
				xx=data.frame(xmeans)  
		}
	}else
	{
		numvar=NULL
		if(!is.null(group.covariates))
			xx=cbind(group.covariates,group=factor(1:nrow(group.covariates)))
		else
			xx=NULL
	}
	chmat=matrix(0,nrow=nocc,ncol=nocc)
	diag(chmat)=1
	ch=apply(chmat,1,paste,collapse=",")
	ch=rep(ch,number.of.groups)
	if(is.null(data$freq))
		data$freq=1
	if(!is.null(xx))
	{
		if(number.of.groups==1)
		{
			data=subset(data,select=c("ch","freq",numvar,names(group.covariates)))
			dummy.data=cbind(data.frame(ch=ch,freq=rep(0,length(ch))),matrix(t(xx),byrow=TRUE,ncol=length(numvar),nrow=length(ch)))
			names(dummy.data)=c("ch","freq",numvar)             
		}else
		{
			data=subset(data,select=c("ch","freq","group",numvar[numvar!="group"],names(group.covariates)))
			xx=subset(xx,select=!names(xx)%in%"group")
			dummy.data=cbind(data.frame(ch=ch,freq=rep(0,length(ch))),group=factor(rep(1:number.of.groups,each=nocc)),
					xx[rep(1:number.of.groups,each=nocc),,drop=FALSE])
			names(dummy.data)=c("ch","freq","group",names(group.covariates),numvar[numvar!="group"])             
		}
	}
	else
	{
		data=subset(data,select=c("ch","freq"))
		dummy.data=data.frame(ch=ch,freq=rep(0,length(ch)))
		names(dummy.data)=c("ch","freq")     
	}
	row.names(dummy.data)=NULL
	return(rbind(data,dummy.data))
}
accumulate_data <- function(data)
{
	x <- data[,names(data)!="freq",drop=FALSE]
	nx <- nrow(x)
	if(is.null(data$freq))data$freq=rep(1,nrow(data))
	pasted.data=apply(x, 1, paste, collapse = "")
	freq=sapply(split(data$freq, pasted.data),sum)
	x=unique(x[order(pasted.data),,drop=FALSE])
	x$freq=freq
	message(nx, " capture histories collapsed into ", nrow(x), "\n", appendLF=FALSE)
	return(x)	
}
