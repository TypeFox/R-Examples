
run.mar<-function(
	data,
	variables=NULL,
	restrictions=NULL,
	search=c("random","exhaustive","fwdstep","exhaustive.true"),
	boot=500,
	ntop=10,
	export=FALSE
	){

search<-search[1]
if(is.na(match(search,c("exhaustive","fwdstep","random","exhaustive.true")))){
stop("The 'search' argument must be equal to \"random\", \"fwdstep\", \"exhaustive\", or \"exhaustive.true\"")
}

run.mar.env<-environment()

if(class(variables)=="MAR") variables<-variables$variables.selected
if(class(restrictions)=="MAR") restrictions<-restrictions$restrictions.set
if(is.null(variables)) restrictions<-NULL

#====================================================================================
# SELECT AND ASSIGN VARIABLES:
#====================================================================================

# Get variable names
heads<-names(data)

if(is.null(variables)){

# Call function for variate/covariate assignment gui
#source("fun_set_variables.R",local=T)

cat(" \n >>>>>> USE BUTTONS TO SELECT EACH VARIATE AND COVARIATE TO INCLUDE IN MODEL <<<<<<
\n      Clicking on each button toggles between: 
            'not included'
            'variate'         (i.e., species or taxa)
            'covariate'       (i.e., exogenous drivers)\n \n \n")
flush.console()

# track holds a sequence of numbers as long as the number of columns
#	in the loaded dataset indicating how each column was assigned
fun.set.variables(heads)->track } else {
variables->track }

if(is.null(track)) {
	on.exit(options(show.error.messages=T))
	options(show.error.messages=F)
	stop(call.=F)}

# Stop if track does not have a value for each column of 'data'
if(length(track)!=length(heads)){
stop("The 'variables' argument must be a numeric vector as long as the number of columns in 'data'.")}
	
# Set variate and covariate matrices and store variable names
year<-data[,1] # dummy variable coding for blocks of data to be considered continuous
time<-data[,2] # time step number

rawvar<-data[,which(track==1),drop=F]  # Variate time-series

# Stop if <1 variates selected
if(ncol(rawvar)<1) {
	variable.error<-"Must select at least 1 variate to include..."
	stop(variable.error)}

rawcovar<-data[,which(track==2),drop=F]  # Covariate time-series

if((ncol(rawvar)+ncol(rawcovar))<2) {
	variable.error<-"Must select at least 1 variate plus 1 covariate or multiple variates to include in model..."
	stop(variable.error)}

namesvar<-heads[which(track==1)]  # save variate names
namescovar<-heads[which(track==2)]  # save covariate names
colnames(rawvar)<-NULL; colnames(rawcovar)<-NULL


#====================================================================================
# SET INTERACTION RESTRICTIONS:
#====================================================================================

if(is.null(restrictions)){

# Call function for variate-variate (B-matrix) and variate-covariate (C-matrix) 
#	interaction restriction gui

cat(" \n >>>>>>>>>>>> USE BUTTONS TO SET RESTRICTIONS ON MODELED INTERACTIONS <<<<<<<<<<<<
\n      Clicking on each button toggles between: 
            0.5     direct interaction is possible (may be included in model)
             0      direct interaction is unlikely/implausible (not included)
             1      direct interaction is probable (will be included in model)\n \n \n")

flush.console()

# Get selected B & C matrix interaction restricitons

fun.restrict.BC(namesvar,namescovar)->indexBCGlobal } else {
restrictions->indexBCGlobal }

{ # 'restrictions' argument ERRORS and WARNINGS
# ERROR: Stop if indexBCGlobal isn't a matrix (wrong restrictions argument class)
if(class(indexBCGlobal)!="matrix"){
stop("The 'restrictions' argument must be of class \"matrix\" or \"MAR\".  See ?run.mar")
}
# ERROR: Stop if indexBCGlobal dimensions do not match number of selected variables
if(sum((dim(indexBCGlobal)-c(length(namesvar),(length(namesvar)+length(namescovar)))))!=0){
stop(paste("The 'restrictions' matrix must have",length(namesvar),"rows and",
	length(namesvar)+length(namescovar),"columns"))
}
# ERROR: Stop if all interactions excluded in restrictions argument
if(sum(indexBCGlobal)==0) stop("All interactions excluded from model!")
# WARNING: Change search type if there are only 1 or 2 models to search through
if(length(which(indexBCGlobal==.5))<2&search!="exhaustive.true"){
search<-"exhaustive.true"
warning("Two or less variable combinations possible.  Model search type \"exhaustive.true\" was used.")
}
}

indexBGlobal<-indexBCGlobal[,1:length(namesvar),drop=F]
if(ncol(indexBGlobal)!=ncol(indexBCGlobal)){
indexCGlobal<-indexBCGlobal[,(length(namesvar)+1):ncol(indexBCGlobal),drop=F]
} else indexCGlobal<-NULL


#====================================================================================
# PREPARE DATA:
#====================================================================================

# Get sizes of matrices
nr<-nrow(rawvar); nc<-ncol(rawvar)
nrc<-nrow(rawcovar); ncc<-ncol(rawcovar)

# Lag the data into separate matrices at time t and time t-1
s1<-nrow(rawvar)
XX<-cbind(rawvar[1:s1-1,],rawvar[2:s1,])

# Get covariates for time t-1
CC<-rawcovar[1:s1-1,,drop=F]

# Throw out non-overlapping time-steps (non-matching 'year' values)
XX<-XX[which(year[1:s1-1]==year[2:s1]),]
lagstate<-as.matrix(XX[,1:nc])					# t-1 variate values
statevar<-as.matrix(XX[,(nc+1):(nc*2)])				# t variate values
covariate<-as.matrix(CC[which(year[1:s1-1]==year[2:s1]),])	# t-1 covariate values

# Initialize variates
Q<-length(statevar[,1])		# length of the time-series
P<-ncol(statevar)			# number of variates in model
R<-ncol(covariate)		# number of covariates in model


#====================================================================================
# INITIALIZE BEST-FIT MODEL SEARCH (lowest AIC):
#====================================================================================

time.start<-Sys.time()		# start timer

# Call the best-fit model search loops
cat(" \n \nsearching for best-fit model...\n")
flush.console()

if(search=="random") search.fun<-bestfit.search.random
if(search=="fwdstep") search.fun<-bestfit.search.fwdstep
if(search=="exhaustive") search.fun<-bestfit.search.exhaustive
if(search=="exhaustive.true") search.fun<-bestfit.search.exhaustive.true

search.fun(statevar,lagstate,covariate,P,R,Q,indexBGlobal,indexCGlobal,ntop)->bestfit.model

#    WAIT...

if(is.null(bestfit.model)) {
	on.exit(options(show.error.messages=T))
	options(show.error.messages=F)
	stop(call.=F)}

bestGlobalB<-bestfit.model$bestGlobalB
bestGlobalC<-bestfit.model$bestGlobalC
all.models<-bestfit.model$all.models

time.end<-Sys.time()		# end timer

elaps<-as.numeric(time.end)-as.numeric(time.start)
mins<-floor(elaps/60)
secs<-format(elaps-60*mins,digits=1)

cat("       ...BEST-FIT MODEL SELECTED
	           ( search time: ",paste(mins,"minutes",secs,"seconds"),")\n \n")


#====================================================================================
# GET TOP BEST-FIT MODELS IF SEARCH TYPE IS 'RANDOM':
#====================================================================================

if(search!="fwdstep"&ntop>1) {
cat(paste(" \n \nidentifying",ntop,"lowest AIC models...\n"))
flush.console()
top.bestfit(ntop,all.models,
	statevar,lagstate,covariate,P,R,Q,namesvar,namescovar)->top.bestfit.models
cat("       ...TOP MODELS RETAINED \n \n")
}

#====================================================================================
# RECALCULATE THE LEAST-SQUARES ESTIMATES FOR THE SELECTED BEST-FIT MODEL:
#====================================================================================

# Run the calculations
bestfit.lstsqr(statevar,lagstate,covariate,P,R,Q,bestGlobalB,bestGlobalC,
						namesvar,namescovar)->best.lstsqr

#====================================================================================
# CALCULATE STABILITY INDICATORS:
#====================================================================================

stab.stats(B=best.lstsqr$B,
			sigma=best.lstsqr$process.errors$covariance,
			V_inf=best.lstsqr$stationary.distribution$covariance)->best.stab


#====================================================================================
# BOOTSTRAPPING:
#====================================================================================

if(boot){

cat(" \n \nbootstrapping best-fit model...\n")
flush.console()
bestfit.bootstrap(best.lstsqr,boot,year,rawvar,rawcovar,s1,P,R,Q,covariate)->bootstrap.results
cat("       ...BOOTSTRAPPING COMPLETE\n \n \n")

# RECALCULATE THE LEAST-SQUARES ESTIMATES FOR THE BOOTSTRAPPED MODEL:

bestfit.lstsqr(statevar,lagstate,covariate,P,R,Q,
			bestGlobalB=bootstrap.results$bootB,
			bestGlobalC=bootstrap.results$bootC,
			namesvar,namescovar)->boot.lstsqr

# CALCULATE STABILITY INDICATORS FOR THE BOOTSTRAPPED MODEL:

stab.stats(B=boot.lstsqr$B,
			sigma=boot.lstsqr$process.errors$covariance,
			V_inf=boot.lstsqr$stationary.distribution$covariance)->boot.stab

}


#====================================================================================
# *** RESULTS ***
#====================================================================================

#cat(paste(rep("=",80),collapse=""),"\n \n \n")
cat(paste(rep("\u2550",80),collapse=""),"\n \n \n")

names(track)<-heads
restrictions<-indexBCGlobal
colnames(restrictions)<-c(namesvar,namescovar)
rownames(restrictions)<-namesvar


MAR.results<-c(
list(
variables.selected	=	track,
restrictions.set		=	restrictions,
search.type		=	search,
search.time.s		=	elaps),
list(bestfit=c(best.lstsqr,stability=list(best.stab)))
)

if(boot){
MAR.results<-c(MAR.results,
	list(bootstrap=c(boot.lstsqr,
		stability=list(boot.stab),
		limits=list(bootstrap.results))
		)
	)
} else MAR.results<-c(MAR.results,list(bootstrap=NULL))

if(search!="fwdstep"&is.numeric(ntop)){
MAR.results<-c(MAR.results,
list(top.bestfit=top.bestfit.models))
}

if(export!=T&export!=F&class(export)!="character") {
	export<-F
	warning("MAR results not exported.  'export' argument must be a character string or logical")
	}

if(export!=F) export.MAR(MAR.results,export)


class(MAR.results)<-"MAR"
MAR.results


}


