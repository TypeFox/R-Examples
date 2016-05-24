#' Simulates data from Hidden Markov Model 
#' 
#' Creates a set of data from a specified HMM model for capture-recapture data.
#' 
#' The specification for the simulation includes a set of data with at least 2 unique ch and freq value to 
#' specify the number of ch values to simulate that start at the specified occasion. For example, 
#' 1000 50
#' 0100 50
#' 0010 50
#' would simulate 150 capture histories with 50 starting at each of occasions 1 2 and 3. The data can also 
#' contain other fields used to generate the model probabilities and each row can have freq=1 to use individual
#' covariates. Either a dataframe (data) is provided and it is processed and the design data list are created
#' or the processed dataframe and design data list are provided. Formula for the model parameters for generating 
#' the data are provided in model.parameters and parameter values are provided in initial.  
#' 
#' @param data Either the raw data which is a dataframe with at least one column named ch (a character field containing the capture history) or a processed dataframe
#' @param ddl Design data list which contains a list element for each parameter type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model 
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional list (by parameter type) of initial values for beta parameters (e.g., initial=list(Phi=0.3,p=-2)
#' @param groups Vector of names of factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce computation
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @export
#' @return dataframe with simulated data
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @examples 
#' # simulate phi(.) p(.) with 1000 Females and 100 males, 3 occasions all released on first occasion
#' df=simHMM(data.frame(ch=c("100","110"),sex=factor(c("F","M")),freq=c(1000,100),
#'    stringsAsFactors=FALSE))
#' df=simHMM(data.frame(ch=rep("100",100),u=rnorm(100,0,1),freq=rep(1,100),
#'   stringsAsFactors=FALSE),
#'   model.parameters=list(Phi=list(formula=~u),p=list(formula=~time)),
#'    initial=list(Phi=c(1,1),p=c(0,1)))
#' df=simHMM(data.frame(ch=c("1000","0100","0010"),freq=rep(50,3),stringsAsFactors=FALSE),
#'   model.parameters=list(Phi=list(formula=~1),p=list(formula=~time)),
#'     initial=list(Phi=c(1),p=c(0,1,2)))
simHMM=function(data,ddl=NULL,begin.time=1,model="hmmCJS",title="",model.parameters=list(),
		design.parameters=list(),initial=NULL,groups=NULL,time.intervals=NULL,accumulate=TRUE,strata.labels=NULL)
{ 
	# call fitHMM to use its code to setup quantities but not fitting model
	setup=crm(data=data,ddl=ddl,begin.time=begin.time,model=model,title=title,model.parameters=model.parameters,
			design.parameters=design.parameters,initial=initial,groups=groups,time.intervals=time.intervals,
			accumulate=accumulate,run=FALSE,strata.labels=strata.labels)
	if(nrow(setup$data$data)==1)stop(" Use at least 2 capture histories; unique ch if accumulate=T")
	parlist=setup$results$par
	T=setup$data$nocc
	m=setup$data$m
	ch=NULL
	df2=NULL
	# compute real parameters
	pars=list()
	for(parname in names(setup$model.parameters))
		pars[[parname]]=do.call("rbind",split(reals(ddl=setup$ddl[[parname]],dml=setup$dml[[parname]],parameters=setup$model.parameters[[parname]],
							parlist=parlist[[parname]]),setup$ddl[[parname]]$id))
	# compute arrays of observation and transition matrices using parameter values
	dmat=setup$data$fct_dmat(pars,m,setup$data$start[,2],T)
	gamma=setup$data$fct_gamma(pars,m,setup$data$start[,2],T)
    delta=setup$data$fct_delta(pars,m,setup$data$start[,2],T,setup$data$start)
	# loop over each capture history
	for (id in as.numeric(setup$data$data$id))
	{
		# set up state with freq rows
		history=matrix(0,nrow=setup$data$data$freq[id],ncol=T)
		state=matrix(0,nrow=setup$data$data$freq[id],ncol=T)
		# create initial state and encounter history value
		state[,setup$data$start[id,2]]=apply(rmultinom(setup$data$data$freq[id],1,delta[id,]),
				                                2,function(x)which(x==1))
		for(k in 1:m)
		{
			instate=sum(state[,setup$data$start[id,2]]==k)
			if(instate>0)
			{
				rmult=rmultinom(instate,1,dmat[id,setup$data$start[id,2],,k])
#				rmult=rmultinom(instate,1,dmat[id,1,,k])
				history[state[,setup$data$start[id,2]]==k,setup$data$start[id,2]]= 
						setup$data$ObsLevels[apply(rmult,2,function(x) which(x==1))]
			}
		}
		# loop over each remaining occasion after the initial occasion 
		for(j in (setup$data$start[id,2]+1):T)
		{
			for(k in 1:m)
			{
				instate=sum(state[,j-1]==k)
				if(instate>0)
				{
					rmult=rmultinom(instate,1,gamma[id,j-1,k,])
					state[state[,j-1]==k,j]= apply(rmult,2,function(x) which(x==1))
				}
			} 
			# use dmat to create observed sequence
			for(k in 1:m)
			{
				instate=sum(state[,j]==k)
				if(instate>0)
				{
					rmult=rmultinom(instate,1,dmat[id,j,,k])
					history[state[,j]==k,j]= setup$data$ObsLevels[apply(rmult,2,function(x) which(x==1))]
				}
			}
		}
	    ch=c(ch,apply(history,1,paste,collapse=","))
		cols2xclude=-which(names(setup$data$data)%in%c("ch","freq","id"))
		if(is.null(df2))
			df2=setup$data$data[rep(id,setup$data$data$freq[id]),cols2xclude,drop=FALSE]
		else
			df2=rbind(df2,setup$data$data[rep(id,setup$data$data$freq[id]),cols2xclude,drop=FALSE])		
	}   
	df=data.frame(ch=ch,stringsAsFactors=FALSE)
	if(nrow(df2)==0)
	   return(df)
    else
	   return(cbind(df,df2))
}

