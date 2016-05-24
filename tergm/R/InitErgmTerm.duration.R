#  File R/InitErgmTerm.duration.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################


InitErgmTerm.edges.ageinterval<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("dissolution","target"))) stop("Term edges.ageinterval can only be used in a dissolution model or as a target statistic.")
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("from","to"),
                      vartypes = c("numeric","numeric"),
                      defaultvalues = list(NULL, Inf),
                      required = c(TRUE, FALSE))
     
  from<-a$from
  to<-a$to

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term edges.ageinterval must have arguments either of the same length, or one of them must have length 1.")
else if(any(from>=to)) stop("Term edges.ageinterval must have from<to.")

  if(any(from<1)) stop("An extant edge cannot have an \"age\" of less than 1.")
  list(name=if(role=="target") "edges_ageinterval_mon" else "edges_ageinterval",
       coef.names = paste("edges.age",from,"to",to,sep=""),
       inputs=c(rbind(from, ifelse(to==Inf, 0, to))),
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm.edge.ages<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term edge.ages can only be used as a target statistic.")
  a <- check.ErgmTerm(nw, arglist)
  
  list(name="edge_ages_mon",
       coef.names = "edge.ages",
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm.mean.age<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term mean.age can only be used as a target statistic.")
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("emptyval","log"),
                      vartypes = c("numeric","logical"),
                      defaultvalues = list(0,FALSE),
                      required = c(FALSE,FALSE))
  
  
  list(name="mean_age_mon",
       coef.names = if(a$log) "mean.log.age" else "mean.age",
       inputs = c(a$emptyval,if(a$log) 1 else 0),
       emptynwstats = a$emptyval,
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm.edgecov.mean.age<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term mean.age can only be used as a target statistic.")

  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname", "emptyval", "log"),
                      vartypes = c("matrix,network,character", "character", "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  
  ### Check the network and arguments to make sure they are appropriate.
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("mean",if(a$log) ".log" else "", ".age", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("mean",if(a$log) ".log" else "", ".age", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(a$emptyval, if(a$log) 1 else 0, NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 3
  list(name="edgecov_mean_age_mon", coef.names = cn, inputs = inputs, duration=TRUE, dependence=FALSE, emptynwstats = a$emptyval)
}

InitErgmTerm.edgecov.ages<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term edgecov.ages can only be used as a target statistic.")

  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network,character", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  
  ### Check the network and arguments to make sure they are appropriate.
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov.ages", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov.ages", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov_ages_mon", coef.names = cn, inputs = inputs, duration=TRUE, dependence=FALSE)
}

################################################################################
InitErgmTerm.degree.mean.age<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term degree.mean.age can only be used as a target statistic.")
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("d", "byarg", "emptyval", "log"),
                      vartypes = c("numeric", "character", "numeric", "logical"),
                      defaultvalues = list(NULL, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$byarg
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degree.mean.age")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree.mean.age() has only one value", call.=FALSE)
  }
  if(!is.null(byarg)) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      stop("Age of ties incident on nodes with degree 0 is meaningless.")
    }
  }
  
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste("degree",d,".mean",if(a$log) ".log" else "", ".age",sep="")
    name <- "degree_mean_age_mon"
    inputs <- c(d)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]],".mean",if(a$log) ".log" else "", ".age", sep="")
    name <- "degree_by_attr_mean_age_mon"
    inputs <- c(as.vector(du), nodecov)
  }
  list(name=name,coef.names=coef.names, inputs=c(a$emptyval, if(a$log) 1 else 0, inputs),
       emptynwstats=rep(a$emptyval,length(coef.names)), duration=TRUE, dependence=TRUE)
}

################################################################################
InitErgmTerm.degrange.mean.age<-function(nw, arglist, role=NULL, ...) {
  if(!is.null(role) && !any(role %in% c("target"))) stop("Term degrange.mean.age can only be used as a target statistic.")
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("from", "to", "byarg", "emptyval", "log"),
                      vartypes = c("numeric", "numeric", "character", "numeric", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, 0, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$byarg
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term edges.ageinterval must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange.mean.age must have from<to.")
  
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degrange.mean.age")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degrange.mean.age() has only one value", call.=FALSE)
  }
  if(!is.null(byarg)) {
    # Combine degree and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      stop("Age of ties incident on nodes with degree 0 is meaningless.")
    }
  }
  
  if(is.null(byarg)) {
    if(length(from)==0){return(NULL)}
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("deg",from,"+",".mean",if(a$log) ".log" else "", ".age",sep=""),
                         paste("deg",from,"to",to,".mean",if(a$log) ".log" else "", ".age",sep=""))
    name <- "degrange_mean_age_mon"
    inputs <- c(rbind(from,to))
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degrange_by_attr function
    coef.names <- ifelse(du[2,]==network.size(nw)+1,
                         paste("deg", du[1,], "+.", byarg, u[du[3,]],".mean",if(a$log) ".log" else "", ".age", sep=""),
                         paste("deg", du[1,], "to", du[2,], ".", byarg, u[du[3,]],".mean",if(a$log) ".log" else "", ".age", sep=""))
    name <- "degrange_by_attr_mean_age_mon"
    inputs <- c(as.vector(du), nodecov)
  }
  list(name=name,coef.names=coef.names, inputs=c(a$emptyval, if(a$log) 1 else 0, inputs),
       emptynwstats=rep(a$emptyval,length(coef.names)), duration=TRUE, dependence=TRUE)
}


