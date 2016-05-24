#  File R/InitErgmTerm.users.R in package ergm.userterms, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
######################################################################
#
# !! USERS: READ THIS FIRST!!
#
# Each term must have its own InitErgmTerm function. You can either
# add all functions to the bottom of this file or have them in 
# separate files (with the extension .R).
#
# This file contains 
#    - a description of the input and output parameters to each
#      InitErgmTerm function
#    - a description of the input parameters that the C changestats
#      function will receive
#    - sample InitErgmTerm functions. These are identical from 
#      "statnet"'s perspective.
#
######################################################################


#  ------------------------------------------------------------------ 
#   Description of the input and output parameters of the  
#   InitErgmTerm.xxx function, where xxx is the name of your term
#  ------------------------------------------------------------------ 
#
#  INPUTS:
#  Each InitErgmTerm function takes three arguments:
#	  		nw: The network of interest
#      arglist: The list of arguments passed to the term xxx
#         ... : There may be other arguments passed by 
#               ergm.getmodel, so each InitErgmTerm function 
#               must include the ... argument
#  These inputs are automatically supplied by ergm.getmodel.
#
#  OUTPUTS:
#  Each InitErgmTerm function should return a list.  
#     REQUIRED LIST ITEMS:
#          name: This names the C changestats function for term xxx, 
#                but does so by excluding the d_ prefix. The 
#                changestats function is named d_xxxy and 'name' is
#                consequently "xxxy". For example, the b1starmix
#                term has 2 changestats functions based on
#                whether the homophily argument is set. These are
#                d_b1starmix and d_b1starmixhomophily. The 'name' 
#                returned by InitErgmTerm.b1starmix is then one of 
#                "b1starmix" or "b1starmixhomophily" as appropriate.
#    coef.names: Vector of names for the coefficients (parameters)
#                as they will be reported in the output.
#       pkgname: This names the package containing the C changestats
#                function d_[name]. The default is "ergm", which means
#                that if you have code that exists as part of the 
#                (say) "ergm.userterms" package, you MUST specify 
#                pkgname="ergm.userterms"
#
#    OPTIONAL LIST ITEMS:
#        inputs: Vector of (double-precision numeric) inputs that the 
#                changestat function called d_'name' may require.
#                The default is NULL; no inputs are required.  But it
#                MUST be a vector!  Thus, if some of the inputs are,  
#                say, matrices, they must be "flattened" to vectors; if 
#                some are categorical character-valued variables, they
#                must be converted to numbers. Optionally, the inputs 
#                vector may have an attribute named "ParamsBeforeCov",
#                which is the number of input parameters preceding the 
#                covariate vector in 'inputs'.  This is necessary for 
#                compatibility with some of the existing d_xxx changestats 
#                functions in ergm, but is not necessary in general.
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  The default value is TRUE.
#  emptynwstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  For example, the degree0 term 
#                would require 'emptynwstats' since degree0 = number of 
#                nodes for the empty network.
#        params: For curved exponential family model terms only, a list of 
#                (numeric) initial values for the parameters of  
#                curved exponential family model terms. Each item in the  
#                list should be named with the corresponding parameter name 
#                (one or more of these will probably coincide with the 
#                 coef.names).  For example, the gwesp term returns 
#                params=list(gwesp=NULL,gwesp.alpha=alpha), where alpha
#                was specified as an argument to the gwesp term. 
#           map: For curved exponential family model terms only, a function 
#                giving the map from the canonical parameters, theta,
#                associated with the statistics for this term, to eta, 
#                the corresponding curved parameters.  The length of eta 
#                is the same as the length of the 'params' list above.
#                The function takes two arguments:  theta and length(eta).
#      gradient: For curved exponential family model terms only, a function 
#                giving the gradient of the 'map'. If theta has length p 
#                and eta has length q, then gradient should return a
#                p by q matrix. This function takes two arguments:  theta 
#                and length(eta).
#


#  ------------------------------------------------------------------------- 
#   Description of the input parameters to the d_xxxy changestats function, 
#   where xxxy corresponds to the 'name' returned by InitErgmTerm.xxx.
#  -------------------------------------------------------------------------- 
#
#  INPUTS:
#  Each d_xxxy function takes five arguments:
#	    ntoggles: the number of toggles as described in 
#                 "ergm.userterms: A template package"
#          heads: a pointer to the array of the head nodes of the 
#                 proposed edges to be toggled
#          tails: a pointer to the array of the tail nodes of the
#                 proposed edges to be toggled
#            mtp: a pointer to the model, which includes the following:
#                 dstats      : a pointer to the array of changestats,
#                               macro-ed as CHANGE_STAT
#                 nstats      : the length of 'dstats', macro-ed as
#                               N_CHANGE_STATS
#                 inputparams : a pointer to the vector of input 
#                               parameters. This is supplied by the
#                               'inputs' returned by InitErgmTerm.xxx
#                               and is macro-ed as INPUT_PARAM
#                 ninputparams: the length of 'inputparams', macro-ed
#                               as N_INPUT_PARAMS
#            nwp: a pointer to the network.  This includes several 
#                 components and several macros exist for accessing
#                 these. See the changestat.h file for a list of these
#                 components and their macros. 
#  These inputs are automatically supplied to the d_xxxy function by the 
#  network_stats_wrapper function 



#  ------------------------------------------------------------------------- 
#  Sample InitErgmTerm function(s)
#  -------------------------------------------------------------------------- 

#  This InitErgmTerm function is for the mindegree term described
#  in Hunter, Goodreau, and Handcock (2011), "ergm.userterms:  A
#  Template Package for Extending statnet"
InitErgmTerm.nmedges<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="nmedges", coef.names="nmedges", pkgnames="ergm.userterms",dependence=FALSE,
       minval = 0, maxval = network.dyadcount(nw), conflicts.constraints="nmedges")
}

InitErgmTerm.nmtriangle<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "nmtriangle")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to triangle() has only one value", call.=FALSE)
    if (!diff) {
      coef.names <- paste("nmtriangle",attrname,sep=".")
      inputs <- c(nodecov)
    } else {
      coef.names <- paste("nmtriangle",attrname, u, sep=".")
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)
    }
  }else{
    coef.names <- "nmtriangle"
    inputs <- NULL
  }
  list(name="nmtriangle", coef.names=coef.names,pkgnames="ergm.userterms", inputs=inputs, minval=0)
}


InitErgmTerm.nmkstar<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("k", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  k<-a$k;attrname<-a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "nmkstar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
#    Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to kstar() has only one value", call.=FALSE)
  }
  lk<-length(k)
  if(lk==0){return(NULL)}
  if(!is.null(attrname)){
    coef.names <- paste("nmkstar",k,".",attrname,sep="")
    inputs <- c(k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lk
  }else{
    coef.names <- paste("nmkstar",k,sep="")
    inputs <- c(k)
  }
  list(name="nmkstar", coef.names=coef.names, pkgnames="ergm.userterms",inputs=inputs, minval = 0, conflicts.constraints="degreedist")
}

InitErgmTerm.blkedges <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("blk","wb","attrname"),
						vartypes = c("numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL),
						required = c(TRUE,TRUE,FALSE))
	wb <- a$wb;blk <- a$blk;attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blkedges")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blkedges() has only one value", call.=FALSE)
	}
	nblk <- length(blk)*(length(blk)+1)/2
	# if(lk==0){
		# return(NULL)
	# }
	if(!is.null(attrname)){
		coef.names <- paste("blkedges",wb,".",attrname,sep="")
		inputs <- c(length(blk),blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- length(blk)+length(wb)+2
	}else{
		coef.names <- paste("blkedges",wb,sep="")
		inputs <- c(length(blk),blk,length(wb),wb)
	}
	list(name="blkedges",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0,conflicts.constraints = "degreedist")
}

InitErgmTerm.blktriangle <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("blk","wb","attrname"),
						vartypes = c("numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL),
						required = c(TRUE,TRUE,FALSE))
	blk <- a$blk;wb <-a$wb;attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blktriangle")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blktriangle() has only one value", call.=FALSE)
	}
	nblk <- length(blk)
	# if(lk==0){
		# return(NULL)
	# }
	if(!is.null(attrname)){
		coef.names <- paste("blktriangle",wb,".",attrname,sep="")
		inputs <- c(length(blk),blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- nblk+length(wb)+2
	}else{
		coef.names <- paste("blktriangle",wb,sep="")
		inputs <- c(length(blk),blk,length(wb),wb)
	}
	list(name="blktriangle",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0)
}

InitErgmTerm.blkkstar <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("k","blk","wb","attrname"),
						vartypes = c("numeric","numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL,NULL),
						required = c(TRUE,TRUE,TRUE,FALSE))
	k <- a$k; blk <- a$blk; wb<-a$wb; attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blkkstar")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blkkstar() has only one value", call.=FALSE)
	}
	nblk <- length(blk)
	lk <-length(k)
	if(lk==0 || nblk==0||length(wb)==0){
		return(NULL)
	}
	coef.names=NULL
	if(!is.null(attrname)){
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
		#coef.names <- paste("k",k,"blkkstar",1:nblk,".",attrname,sep="")
		inputs <- c(lk,k,nblk,blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- c(lk+length(blk)+length(wb)+3)
		#attr(inputs,"ParamsBeforeblk") <- lk
	}else{
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
	#	coef.names <- paste("k",k,"blkkstar",1:nblk,sep="")
		inputs <- c(lk,k,nblk,blk,length(wb),wb)
	}
	list(name="blkkstar",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0,conflicts.constraints = "degreedist")
}

#the init term function for the second type of block model.
InitErgmTerm.blktriangle2 <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("blk","wb","attrname"),
						vartypes = c("numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL),
						required = c(TRUE,TRUE,FALSE))
	blk <- a$blk;wb <- a$wb;attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blktriangle2")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blktriangle2() has only one value", call.=FALSE)
	}
	nblk <- length(blk)
	# if(lk==0){
		# return(NULL)
	# }
	if(!is.null(attrname)){
		coef.names <- paste("blktriangle",wb,".",attrname,sep="")
		inputs <- c(nblk,blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- nblk+wb+2
	}else{
		coef.names <- paste("blktriangle",wb,sep="")
		inputs <- c(nblk,blk,length(wb),wb)
	}
	list(name="blktriangle2",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0)
}


# InitErgmTerm.blktriangle2 <- function(nw,arglist,...){
	# a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						# varnames = c("blk","attrname"),
						# vartypes = c("numeric","character"),
						# defaultvalues = list(NULL,NULL),
						# required = c(TRUE,FALSE))
	# blk <- a$blk;attrname <- a$attrname
	# if(!is.null(attrname)){
		# nodecov <- get.node.attr(nw, attrname, "blktriangle2")
    	# u<-sort(unique(nodecov))
	    # if(any(is.na(nodecov))){u<-c(u,NA)}
	# #    Recode to numeric if necessary
    	# nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    # if (length(u)==1)
    	  # stop ("Attribute given to blktriangle2() has only one value", call.=FALSE)
	# }
	# nblk <- length(blk)
	# # if(lk==0){
		# # return(NULL)
	# # }
	# if(!is.null(attrname)){
		# coef.names <- paste("blktriangle",1:nblk,".",attrname,sep="")
		# inputs <- c(blk,nodecov)
		# attr(inputs,"ParamsBeforeCov") <- nblk
	# }else{
		# coef.names <- paste("blktriangle",1:nblk,sep="")
		# inputs <- c(blk)
	# }
	# list(name="blktriangle2",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0)
# }


InitErgmTerm.blkkstar2 <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("k","blk","wb","attrname"),
						vartypes = c("numeric","numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL,NULL),
						required = c(TRUE,TRUE,TRUE,FALSE))
	k <- a$k; blk <- a$blk; wb <-a$wb; attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blkkstar2")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blkkstar2() has only one value", call.=FALSE)
	}
	nblk <- length(blk)
	lk <-length(k)
	if(lk==0 || nblk==0||length(wb)==0){
		return(NULL)
	}
	coef.names=NULL
	if(!is.null(attrname)){
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
		#coef.names <- paste("k",k,"blkkstar",1:nblk,".",attrname,sep="")
		inputs <- c(lk,k,nblk,blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- lk+nblk+length(wb)+3
		#attr(inputs,"ParamsBeforeblk") <- lk
	}else{
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
	#	coef.names <- paste("k",k,"blkkstar",1:nblk,sep="")
		inputs <- c(lk,k,nblk,blk,length(wb),wb)
	}
	list(name="blkkstar2",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0,conflicts.constraints = "degreedist")
}


InitErgmTerm.blktriangle3 <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("blk","wb","attrname"),
						vartypes = c("numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL),
						required = c(TRUE,TRUE,FALSE))
	blk <- a$blk;wb <-a$wb; attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blktriangle3")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blktriangle2() has only one value", call.=FALSE)
	}
	nblk <- length(blk)*(length(blk)+1)/2
	# if(lk==0){
		# return(NULL)
	# }
	if(!is.null(attrname)){
		coef.names <- paste("blktriangle",wb,".",attrname,sep="")
		inputs <- c(length(blk),blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- length(blk)+length(wb)+2
	}else{
		coef.names <- paste("blktriangle",wb,sep="")
		inputs <- c(length(blk),blk,length(wb),wb)
	}
	list(name="blktriangle3",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0)
}

InitErgmTerm.blkkstar3 <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("k","blk","wb","attrname"),
						vartypes = c("numeric","numeric","numeric","character"),
						defaultvalues = list(NULL,NULL,NULL,NULL),
						required = c(TRUE,TRUE,TRUE,FALSE))
	k <- a$k; blk <- a$blk; wb <-a$wb; attrname <- a$attrname
	if(!is.null(attrname)){
		nodecov <- get.node.attr(nw, attrname, "blkkstar3")
    	u<-sort(unique(nodecov))
	    if(any(is.na(nodecov))){u<-c(u,NA)}
	#    Recode to numeric if necessary
    	nodecov <- match(nodecov,u,nomatch=length(u)+1)
	    if (length(u)==1)
    	  stop ("Attribute given to blkkstar2() has only one value", call.=FALSE)
	}
	nblk <- length(blk)*(length(blk)-1)/2
	lk <-length(k)
	if(lk==0 || nblk==0 ||length(wb)==0){
		return(NULL)
	}
	coef.names=NULL
	if(!is.null(attrname)){
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
		#coef.names <- paste("k",k,"blkkstar",1:nblk,".",attrname,sep="")
		inputs <- c(lk,k,length(blk),blk,length(wb),wb,nodecov)
		attr(inputs,"ParamsBeforeCov") <- c(lk+length(blk)+length(wb)+3)
		#attr(inputs,"ParamsBeforeblk") <- lk
	}else{
		for(i in 1:length(wb))
		{
			for(j in 1:lk)
			{
				coef.names <- c(coef.names,paste("blk",wb[i],"kstar",k[j],".",attrname,sep=""))
			}
		}
	#	coef.names <- paste("k",k,"blkkstar",1:nblk,sep="")
		inputs <- c(lk,k,length(blk),blk,length(wb),wb)
	}
	list(name="blkkstar3",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0,conflicts.constraints = "degreedist")
}

InitErgmTerm.degseq<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = FALSE)
     coef.names <- paste("node",1:nw$gal$n,"deg",sep="")
  list(name="degseq", coef.names=coef.names, pkgnames="ergm.userterms",minval = 0, conflicts.constraints="degreedist")
}


InitErgmTerm.blkdegseq <- function(nw,arglist,...){
	a <- check.ErgmTerm(nw,arglist,directed=FALSE,
						varnames = c("blk","wb"),
						vartypes = c("numeric","numeric"),
						defaultvalues = list(NULL,NULL),
						required = c(TRUE,TRUE))
	blk <- a$blk;wb <- a$wb
	nblk <- length(blk)
	coef.names=NULL
	for(i in 1:length(wb))
	{
		coef.names <- c(coef.names,paste("blk",wb[i],"degseq",(blk[wb[i]]+1):blk[wb[i]+1],sep=""))
	}
	inputs <- c(nblk,blk,length(wb),wb)
	list(name="blkdegseq",coef.names=coef.names,pkgnames="ergm.userterms",inputs=inputs,minval=0,conflicts.constraints = "degreedist")
}








