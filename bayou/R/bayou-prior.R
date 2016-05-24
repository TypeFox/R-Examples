#' Make a prior function for bayou
#' 
#' This function generates a prior function to be used for bayou according to user specifications.
#'
#' @param tree A tree object of class "phylo"
#' @param dists A list providing the function names of the distribution functions describing the prior distributions of parameters (see details). If no
#' distributions are provided for a parameter, default values are given. Note that the names are provided as text strings, not the functions themselves.
#' @param param A list providing the parameter values of the prior distributions (see details). 
#' @param plot.prior A logical indicating whether the prior distributions should be plotted.
#' @param model One of three specifications of the OU parameterization used. 
#' Takes values \code{"OU"} (alpha & sig2), \code{"QG"} (h2, P, w2, Ne), or \code{"OUrepar"} (halflife,Vy)
#' @param fixed A list of parameters that are to be fixed at provided values. These are removed from calculation of the prior value. 
#' @details Default distributions and parameter values are given as follows:
#' OU: \code{list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm",
#'          "dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),
#'    param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),
#'          "dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
#' QG: \code{list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm",
#'          "dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),
#'    param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),
#'          "dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list()))}
#' OUrepar: \code{list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm",
#'          "dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),
#'    param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),
#'          "dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1)),"dsb"=list())}
#' 
#' \code{dalpha, dsig2, dh2, dP, dw2, dNe, dhalflife},  and \code{dVy} must be positive continuous distributions and provide the parameters used to calculate alpha and sigma^2 of the OU model. 
#' \code{dtheta} must be continuous and describes the prior distribution of the optima. dk is the prior distribution for the number of shifts. For Poisson and conditional Poisson (cdpois) are provided
#' the parameter \code{lambda}, which provides the total number of shifts expected on the tree (not the rate per unit branch length). Otherwise, \code{dk} can take any positive, discrete distribution.  
#' dsb indicates the prior probability of a given set of branches having shifts, and is generally specified by the "dsb" function in the bayou package. See the documentation for dsb for specifying the number
#' of shifts allowed per branch, the probability of a branch having a shift, and specifying constraints on where shifts can occur.\code{"dloc"} indicates the prior probability of the location of a shift within
#' a single branch. Currently, all locations are given uniform density. All distributions are set to return log-transformed probability densities. 
#' 
#' @return returns a prior function of class "priorFn" that calculates the log prior density for a set of parameter values provided in a list with correctly named values.
#' 
#' @export
#' @examples 
#' ## Load data
#' data(chelonia)
#' tree <- chelonia$phy
#' dat <- chelonia$dat
#' 
#' #Create a prior that allows only one shift per branch with equal probability 
#' #across branches
#' prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="dlnorm",
#'            dsb="dsb", dk="cdpois", dtheta="dnorm"), 
#'              param=list(dalpha=list(meanlog=-5, sdlog=2), 
#'                dsig2=list(meanlog=-1, sdlog=5), dk=list(lambda=15, kmax=200), 
#'                  dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(dat), sd=2)))
#'              
#' #Evaluate some parameter sets
#' pars1 <- list(alpha=0.1, sig2=0.1, k=5, ntheta=6, theta=rnorm(6, mean(dat), 2), 
#'                  sb=c(32, 53, 110, 350, 439), loc=rep(0.1, 5), t2=2:6)
#' pars2 <- list(alpha=0.1, sig2=0.1, k=5, ntheta=6, theta=rnorm(6, mean(dat), 2),
#'                  sb=c(43, 43, 432, 20, 448), loc=rep(0.1, 5), t2=2:6)
#' prior(pars1) 
#' prior(pars2) #-Inf because two shifts on one branch
#' 
#' #Create a prior that allows any number of shifts along each branch with probability proportional 
#' #to branch length
#' prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="dlnorm",
#'            dsb="dsb", dk="cdpois", dtheta="dnorm"), 
#'              param=list(dalpha=list(meanlog=-5, sdlog=2), 
#'                dsig2=list(meanlog=-1, sdlog=5), dk=list(lambda=15, kmax=200), 
#'                  dsb=list(bmax=Inf,prob=tree$edge.length), 
#'                    dtheta=list(mean=mean(dat), sd=2)))
#' prior(pars1)
#' prior(pars2)
#' 
#' #Create a prior with fixed regime placement and sigma^2 value
#' prior <- make.prior(tree, dists=list(dalpha="dlnorm", dsig2="fixed", 
#'            dsb="fixed", dk="fixed", dtheta="dnorm", dloc="dunif"), 
#'              param=list(dalpha=list(meanlog=-5, sdlog=2), 
#'                dtheta=list(mean=mean(dat), sd=2)), 
#'                  fixed=list(sig2=1, k=3, ntheta=4, sb=c(447, 396, 29)))
#'                  
#' pars3 <- list(alpha=0.01, theta=rnorm(4, mean(dat), 2), loc=rep(0.1, 4))
#' prior(pars3)
#' 
#' ##Return a list of functions used to calculate prior
#' attributes(prior)$functions
#' 
#' ##Return parameter values used in prior distribution
#' attributes(prior)$parameters

make.prior <- function(tree, dists=list(), param=list(), fixed=list(), plot.prior=TRUE,model="OU"){
  tree <- reorder.phylo(tree, "postorder")
  nH <- max(nodeHeights(tree))
  ntips <- length(tree$tip.label)
  TH <- sum(tree$edge.length)
  default.OU <- list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dalpha"=list(),"dsig2"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  default.QG <- list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dtheta"=list(),"dk"=list(lambda=1,kmax=2*ntips-2),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  default.OUrepar <- list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dk"="cdpois","dtheta"="dnorm","dsb"="dsb","dloc"="dunif"),param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1,kmax=2*ntips-2),"dtheta"=list(),"dloc"=list(min=0,max=1),"dsb"=list(ntips=ntips, bmax=1, prob=1)))
  #default.OUcpp <- list(dists=list("dalpha"="dlnorm","dsig2"="dlnorm","dsig2jump"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dalpha"=NULL,"dsig2"=list(),"dsig2jump"=list(),"dtheta"=list(),"dk"=list(lambda=1),"dloc"=list(min=0,max=TH)))
  #default.QGcpp <- list(dists=list("dh2"="dbeta","dP"="dlnorm","dw2"="dlnorm","dNe"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dh2"=list(shape1=1,shape2=1),"dP"=list(),"dw2"=list(),"dNe"=list(),"dsig2jump"=list(),"dtheta"=list(),"dk"=list(lambda=1),"dloc"=list(min=0,max=TH)))
  #default.OUreparcpp <- list(dists=list("dhalflife"="dlnorm","dVy"="dlnorm","dsig2jump"="dlnorm","dk"="dpois","dtheta"="dnorm","dloc"="dunif"),param=list("dhalflife"=list("meanlog"=0.25,"sdlog"=1.5),"dVy"=list("meanlog"=1,"sdlog"=2),"dk"=list(lambda=1),"dsig2jump"=list(),"dtheta"=list(),"dloc"=list(min=0,max=TH)))
  default <- switch(model,"OU"=default.OU,"QG"=default.QG,"OUrepar"=default.OUrepar)#,"OUcpp"=default.OUcpp,"QGcpp"=default.QGcpp,"OUreparcpp"=default.OUreparcpp)
  notprovided <- setdiff(names(default$dist),names(dists))
  pars.notprovided <- setdiff(names(default$param),names(param))
  dists[notprovided] <- default$dists[notprovided]
  param[pars.notprovided] <- default$param[pars.notprovided]
  if(length(setdiff(names(default$param),names(param)))>0)
    stop("Provided parameters are not in the model")
  if(length(setdiff(names(default$dists),names(dists)))>0)
    stop("Provided parameters are not in the model")
  if(dists$dsb=="dsb") param$dsb$ntips <- ntips
  if(dists$dloc=="dunif" & dists$dsb=="dsb"){
        dists$dloc <- "dloc"
    }
  remove <- which(unlist(dists)=="fixed")
  if(length(remove)>0){
    dists2get <- dists[-remove]
    param2get <- param[-(match(names(remove),names(param)))]
    fixed.param <- gsub('^[a-zA-Z]',"",names(remove)) 
    if("sb" %in% fixed.param & !("k" %in% names(fixed))) fixed$k <- length(fixed$sb)
    if("sb" %in% fixed.param & !("t2" %in% names(fixed)) & length(fixed$sb)>0) fixed$t2 <- 2:(length(fixed$sb)+1)
    missing.fixed <- fixed.param[!(fixed.param %in% names(fixed))]
    if(length(missing.fixed)>0){
      stop(paste("'", paste(missing.fixed, collapse="', '"), "' set as 'fixed' but not provided", sep=""))
    }
  } else {dists2get <- dists; param2get <- param}
  prior.fx <- lapply(dists2get,get)
  param2get <- suppressWarnings(lapply(param2get,function(x){ x$log = TRUE; x}))
  prior.param <- param2get[match(names(prior.fx),names(param2get))]
  prior.fx <- lapply(1:length(prior.param),function(x) .set.defaults(prior.fx[[x]],defaults=prior.param[[x]]))
  names(prior.fx) <- names(prior.param)
  #if(model %in% c("OUcpp","QGcpp","OUreparcpp")){
    #droot <- prior.fx$dtheta
    #prior.fx$dtheta <- function(x){
      #droot(x[1])
    #}
  #}
  par.names <- gsub('^[a-zA-Z]',"",names(dists2get))
  
  if(plot.prior){
    par(mfrow=c(ceiling(length(dists2get)/2),2))
    rfx <- lapply(gsub('^[a-zA-Z]',"r",dists2get),function(x) try(get(x),silent=TRUE))
    rprior.param <- prior.param[1:(length(prior.param))]
    rprior.param <- lapply(rprior.param, function(x) x[-length(x)])
    if(!is.null(dists2get$dsb) & !is.null(dists2get$dloc)){
      if(dists2get$dsb=="dsb" & any(rprior.param$dsb$bmax==1)){rprior.param$dsb$bmax[rprior.param$dsb$bmax==1] <- Inf; rprior.param$dsb$prob <- 1}
    }
    rfx <- lapply(1:length(rprior.param),function(x) try(.set.defaults(rfx[[x]],defaults=rprior.param[[x]]),silent=TRUE))
    plot.names<-par.names[sapply(rfx,class)=="function"]
    rfx <- rfx[sapply(rfx,class)=="function"]
    names(rfx) <- names(rprior.param)
    nsim <-500000
    for(i in 1:length(rfx)){
      if(names(rfx)[i]=="dsb"){
        curve(sapply(x,function(y) prior.fx$dsb(y,log=FALSE)),xlim=c(1,(2*ntips-2)),ylab="Density",main="branches")
      } else {
          x <- rfx[[i]](nsim)
          qq <- quantile(x,c(0.001,0.999))
          plot(density(x),xlim=qq, main=plot.names[i],lwd=2)
          }
      }  
  }
  priorFUN <- function(pars,cache){
      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
      pars.o <- pars[match(par.names,names(pars))]
      pars.o <- pars.o[!is.na(names(pars.o))]
      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
      names(densities) <- par.names
      lnprior <- sum(unlist(densities,F,F))
      return(lnprior)
  }

#  if(type=="emap"){
#    priorFUN <- function(pars,cache,emap){
#      pars$sb <- which(emap$sh==1)
#      pars$loc <- emap$r1[pars$sb]
#      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
#      pars.o <- pars[match(par.names,names(pars))]
#      pars.o <- pars.o[!is.na(names(pars.o))]
#      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
#      names(densities) <- par.names
#      lnprior <- sum(unlist(densities,F,F))
#      return(lnprior)
#    }
#  }
#  if(type=="simmap"){
#    priorFUN <- function(pars,cache){
#      pars$sb <- rep(1:length(cache$edge.length),sapply(cache$maps,length)-1)
#      pars$loc <- sapply(pars$sb,function(x) cache$maps[[x]][-1])
#      pars$t2 <- names(pars$loc)
#      pars$loc <- unname(pars$loc)
#      if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
#      pars.o <- pars[match(par.names,names(pars))]
#      pars.o <- pars.o[!is.na(names(pars.o))]
#      densities <- sapply(1:length(pars.o),function(x) prior.fx[[x]](pars.o[[x]]))
#      names(densities) <- par.names
#      lnprior <- sum(unlist(densities,F,F))
#      return(lnprior)
#    }
#  }
  if(length(remove)>0){
    if("sb" %in% fixed.param){
      prior.param$dsb$bmax <- rep(0, nrow(tree$edge))
      prior.param$dsb$bmax[fixed$sb] <- 1
      prior.param$dsb$prob <- prior.param$dsb$bmax
    }
  }
  attributes(priorFUN) <- list("model"=model,"parnames"=par.names,"distributions"=dists,"parameters"=prior.param,"fixed"=fixed, "functions"=prior.fx)
  class(priorFUN) <- c("priorFn","function")
  return(priorFUN)
}
