validate <- function ( generate.param, generate.param.inputs=NULL, 
	generate.data, generate.data.inputs=NULL, analyze.data, 
	analyze.data.inputs=NULL, n.rep=20,n.batch=NULL,params.batch=NULL, 
	print.reps=FALSE) {

#----------------------------------------------------------------------------
#        Inputs
#
# generate.param : function for generating parameters from prior distribution
#                  should output a vector of parameters
#                  function should look like generate.param <- function() {...}
#                  or generate.param <- function(generate.param.inputs) {...}
# generate.data  : function for generating data given parameters
#                  should take as input the output from generate.param
#                  should output data matrix to be analyzed
#                  function should look like 
#		   generate.data <- function(theta.true) {...} or
#                  generate.data <- function(theta.true, generate.data.inputs)
#		   {...}
# analyze.data   : function for generating sample from posterior distribution 
#                  should take as input the output from generate.data and 
#		   generate.param 
#                  should output a matrix of parameters, each row is one 
#		   parameter draw
#                  function should look like 
#		   analyze.data <- function(data.rep,theta.true) {...} or 
#                  analyze.data <- function(data.rep,theta.true, 
#		   analyze.data.inputs) {...}
#                !! This is the software being tested !!
#
#----------------------------------------------------------------------------
#
#        Optional Inputs
#
# generate.param.inputs : inputs to the function generate.param
# generate.data.inputs  : inputs to the function generate.data 
#			  (in addition to theta.true)
# analyze.data.inputs   : inputs to the function analyze.data 
#			  (in addition to data.rep and theta.true)
# n.batch               : lengths of parameter batches.  Must sum to n.param 
#			  and correspond to the order of parameters as they 
#			  are output from analyze.data.  For example, if 
#			  there are 5 total parameters with the first two in 
#			  one batch and the next three in another batch, 
#			  use n.batch=c(2,3)
# params.batch		: names of parameter batches, used as the y axis in 
#			  the output plot.  Must have length equal to the 
#			  number of batches.  Can consist of text (e.g., 
#			  params.batch=c("alpha","beta")) or an expression 
#			  (e.g., params.batch=expression(alpha,beta)).  
#			  Not used if n.batch not provided.
# n.rep                 : number of replications to be performed, default is 20
# print.reps            : indcator of whether or not to print the replication 
#			  number, default is FALSE
#
#-----------------------------------------------------------------------------
#
#        Variables
#
# n.param           : dimension of parameter, total number of parameters, 
#		      length of theta
# theta.true        : vector of true parameters, used to generate data and 
#		      calculate posterior quantiles.  length equals n.param 
# data.rep          : matrix of data generated from theta.true
# theta.draws       : sample from posterior distribution
#                     matrix with number of columns equals n.param
# n.draws           : number of rows in theta.draws, size of posterior 
#		      sample generated
# quantiles.theta   : matrix of posterior quantiles of each component of 
#		      theta.true, number of rows equals n.rep, number of 
#		      columns equals n.param
# p.vals	    : p-values based on testing for uniformity of the 
#		      quantiles of the components of each scalar component 
#		      of theta
# p.batch	    : p-values for parameter "batches", e.g., means of 
#	              person-level parameters.  These are the p-values used 
#                     in the Bonferroni correction for multiple comparisons
# adj.min.p         : the smallest of the batched p-values, with Bonferroni
#                     correction applied
# z.stats	    : z statistics corresponding to p-vals and p.batch, 
#		      for plotting

#----------------------------------------------------------------------------
#            Initial Setup

##see how long the parameter vector is
if(is.null(generate.param.inputs)) n.param <- length(generate.param()) else 
	n.param <-length(generate.param(generate.param.inputs))

##create variables for batched parameters
if(!is.null(n.batch)) {
	num.batches<-length(n.batch)
  ##first, make sure batch parameter lengths match with parameter vector length
	if(sum(n.batch) != n.param){ 
	print ("Error:  Lengths of parameter batches don't sum to length of parameter vector")
	return()}
	if(!is.null(params.batch)){
	if(length(params.batch) != length(n.batch)){
		print("Error:  Must have same number of named parameters as batches")
		return()}
	}
  ##indices to take means over
	batch.ind <- rep(0,(num.batches+1))
	for(i in 1:num.batches) batch.ind[(i+1)] <- batch.ind[i] + n.batch[i]
  ##axis values for plot
	plot.batch <- rep(1,n.batch[1])
	for(i in 2:num.batches) plot.batch<-c(plot.batch,rep(i,n.batch[i]))
	}

##initialize quantiles
##if parameters are in batches, make extra columns in quantile matrix for 
##batch means
if(!is.null(n.batch)) 
	quantile.theta <- matrix(0,nrow=n.rep,ncol=(n.param+num.batches)) else
	quantile.theta <- matrix(0,nrow=n.rep,ncol=n.param)

#---------------------------------------------------------------------------
#            Validation Loop

for(reps in 1:n.rep){

	##First generate parameters
	##theta.true is vector of true parameters, with length n.param
	if(is.null(generate.param.inputs)) theta.true <- generate.param() else 
		theta.true<-generate.param(generate.param.inputs)

	##Next generate data
        if(is.null(generate.data.inputs)) 
		data.rep <- generate.data(theta.true) else
		data.rep <- generate.data(theta.true, generate.data.inputs)

	##Now analyze data
	if(is.null(analyze.data.inputs)) 
		theta.draws <- analyze.data(data.rep,theta.true) else
		theta.draws <- analyze.data(data.rep,theta.true, analyze.data.inputs)

	if ( is.matrix(theta.draws) == FALSE )  
		theta.draws<-as.matrix(theta.draws)

	##Make sure that theta.true and theta.draws have same number 
	##of parameters
	if( length(theta.true) != ncol(theta.draws) ){
	   print("Error: Generate.param and analyze.data must output same number of parameters")
	   return()}

	##If parameters are batched, add means of batched parameters to 
	##output matrix
	if(!is.null(n.batch)){
	  for(i in 1:num.batches) {
	    if(n.batch[i]>1){
	    theta.draws <- cbind(theta.draws,
	      apply(as.matrix(theta.draws[,(batch.ind[i]+1):batch.ind[(i+1)]]),
	      1,mean))
	    theta.true <- c(theta.true,
	      mean(theta.true[(batch.ind[i]+1):batch.ind[(i+1)]]))
	    } else {
	    theta.draws <- cbind(theta.draws,theta.draws[,(batch.ind[i]+1)])
	    theta.true <- c(theta.true, theta.true[(batch.ind[i]+1)])
	    }
	  }
	}

	##update matrix of posterior quantiles
	theta.draws <- rbind(theta.true, theta.draws)
	quantile.theta[reps, ] <- apply( theta.draws, 2, quant )

	if(print.reps==TRUE) print(reps)
}


#-----------------------------------------------------------------------------
#            Analyze simulation output

##calculate z.theta statistics and p-values
quantile.trans <- (apply(quantile.theta, 2, qnorm))^2
q.trans <- apply(quantile.trans,2,sum)
p.vals <- pchisq(q.trans,df=n.rep,lower.tail=FALSE)
z.stats <- abs( qnorm(p.vals) )

lower.lim<-min( min(z.stats-1), -3.5)
upper.lim<-max( max(z.stats+1), 3.5)

##Bonferonni correction
if (is.null(n.batch)) {
    adj.min.p <- n.param*min(p.vals)
}
else {
    z.batch <- z.stats[(n.param + 1):length(p.vals)]
    p.batch <- p.vals[(n.param + 1):length(p.vals)]
    adj.min.p <- num.batches*min(p.batch)
}

##Print output
print(paste("Smallest Bonferonni-adjusted p-value: ", round(adj.min.p, 3)))

##plot
if(is.null(n.batch)){
	plot(z.stats, rep(1,n.param), xlim=c(0,upper.lim),xlab="",
		ylab="",main=expression("Absolute "*z[theta]*" Statistics"),
		axes=F)
	axis(1,line=.1)} else {
	##first plot parameters that are NOT batch means
	rows.one<-c(1:num.batches)[n.batch==1]
	plot(z.stats[1:n.param][plot.batch%in%rows.one==FALSE],
		plot.batch[plot.batch%in%rows.one==FALSE],
		xlim=c(0,upper.lim),ylim=c(1,num.batches),xlab="",
		ylab="",main=expression("Absolute "*z[theta]*" Statistics"),
		axes=F,pch=176)
	##now add batch means
	points(z.batch,c(1:num.batches),pch=20)
	axis(1,line=.1)
	if(!is.null(params.batch)) 	
		axis(2,at=c(1:max(plot.batch)),tick=FALSE,labels=params.batch,
		las=1,line=-1)}
abline(v=0)

#----------------------------------------------------------------------------

##return z statistics and p-value
if(is.null(n.batch)){
	return(list(p.vals = p.vals, adj.min.p=adj.min.p))} else {
	
	if(length(z.batch)==n.param)
	return(list(p.batch = p.batch, adj.min.p=adj.min.p)) else 
	return(list(p.vals = p.vals[1:n.param], p.batch = p.batch, 
            adj.min.p=adj.min.p))
	}

}


