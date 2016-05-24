#ALG: this code is based on the code 'rpart.control' but
# modifies/extends it with the additional procedures and their
# parameters.

#SCCS @(#)rpart.control.s	1.10 07/05/01
# alg 2/1/2012: added default for penalty term in interpretable trees
# alg 2/11/2012: added inter_param=c(0,0)
# alg 4/19/2012: added improve scaling control
#
itree.control <-
  function(minsplit=20, minbucket= round(minsplit/3),
	   maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
	   surrogatestyle =0, maxdepth=30,impscale=3,
	   interp_param1=0,interp_param2=0,cp,... ) {

	extraArgs = list(...)

	if (maxcompete<0) {
	    warning("The value of maxcompete supplied is <0; the value 0 was used instead")
	    maxcompete <-0
	    }
	if (any(xval<0)) {
	    warning("The value of xval supplied is <0; the value 0 was used instead")
	    xval <-0
	    }
	if (maxdepth > 30) stop("Maximum depth is 30")
	if (maxdepth < 1)  stop("Maximum depth must be at least 1")

	if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket*3

	if((usesurrogate < 0) || (usesurrogate > 2)) {
	    warning("The value of usesurrogate supplied was out of range," ,
		    "the default value of 2 is used instead.")
	    usesurrogate <- 2
	    }
	if((surrogatestyle < 0) || (surrogatestyle > 1)) {
	    warning("The value of surrogatestyle supplied was out of range,",
		    "the default value of 0 is used instead.")
	    surrogatestyle <- 0
	    }

	#ALG 4/19
	if(impscale!=1 && impscale!=2 && impscale !=3){
		warning("The value of impscale was not 1 (parent) or 2 (root), ",
				"the default for this combination of splitting objective and penalty will be used.")
		impscale <- 3
	}

	#ALG 4/27. if cp not passed, set to default for the method in rpart.s code.
	if(missing(cp)){
		cp <- NULL
	}

	# Because xval can be of length either 1 or n, and the C code
	#   refers to parameters by number, i.e., "opt[5]" in rpart.c,
	#   the xval parameter should always be last on the list.
	#
	# alg 2/11/2012: added interp_param1 and 2
	# alg 4/19/2012: added impscale
	list(minsplit=minsplit, minbucket=minbucket, cp=cp,
	     maxcompete=maxcompete, maxsurrogate=maxsurrogate,
	     usesurrogate=usesurrogate,
	     surrogatestyle=surrogatestyle, maxdepth=maxdepth,
	     impscale=impscale,
	     interp_param1=interp_param1,
	     interp_param2=interp_param2, xval=xval)
	}
