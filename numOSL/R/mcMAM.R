#####
mcMAM<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    UseMethod("mcMAM")
} #
### 2014.09.19.
mcMAM.default<-
function(EDdata, ncomp=-1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    ### Stop if not.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,2L,drop=TRUE]>0),
              length(ncomp)==1L, ncomp %in% c(-1L,-2L),
              length(addsigma)==1L, is.numeric(addsigma),
              length(iflog)==1L, is.logical(iflog),
              length(nsim)==1L, is.numeric(nsim), nsim>=100L, nsim<=2e5,
              is.list(inis), is.list(control.args),
              all(names(inis) %in% c("p","gamma","mu","sigma")), 
              all(names(control.args) %in% c("w","m","nstart")))
    ###
    ed1<-as.numeric(EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(EDdata[,2L,drop=TRUE])
    nED<-nrow(EDdata)
    ###
    if (iflog==TRUE && any(ed1<=0)) {
        stop("Error: minus(zero) ED canot be logged!")
    } # end if
    ### Bounds of dose population.
    rangeED<-range(ed1)
    if(all(ed1>0))  {
        lowerGamma<-rangeED[1L]*0.999
        upperGamma<-rangeED[2L]*1.001
    } else if (all(ed1<=0))  {
        lowerGamma<-rangeED[1L]*1.001
        upperGamma<-rangeED[2L]*0.999
    } else {
        lowerGamma<-rangeED[1L]*1.001
        upperGamma<-rangeED[2L]*1.001
    } # end if
    ###
    ### Default inis.
    args.inis<-list("p"=0.5,"gamma"=quantile(ed1,probs=0.25,names=FALSE),
                    "mu"=quantile(ed1,probs=0.5,names=FALSE),"sigma"=0.5)
    args.inis[names(inis)]<-inis
    stopifnot(args.inis[["p"]]>0, args.inis[["p"]]<1,
              args.inis[["gamma"]]>lowerGamma, 
              args.inis[["gamma"]]<upperGamma,
              args.inis[["mu"]]>lowerGamma, 
              args.inis[["mu"]]<upperGamma,
              args.inis[["sigma"]]>0,
              args.inis[["sigma"]]<ifelse(iflog==TRUE,5,var(ed1)))
    ###
    ### Default arguments of slice sampling.
    args.control<-list(w=1, m=-100, nstart=1L)
    args.control[names(control.args)]<-control.args
    stopifnot(args.control[["w"]]>=1e-2,
              args.control[["w"]]<=1e3,
              args.control[["m"]]<=1e9,
              args.control[["nstart"]]>=1L,
              args.control[["nstart"]]<=1000L)
    ###
    w<-args.control$w
    m<-args.control$m
    nstart<-args.control$nstart
    inis<-if (ncomp==-1L) {
        c(args.inis[["p"]],
          args.inis[["gamma"]],
          args.inis[["sigma"]])
    } else {
        c(args.inis[["p"]],
          args.inis[["gamma"]],
          args.inis[["mu"]],
          args.inis[["sigma"]])
    } # end if
    iflag<-0
    chains<-matrix(0,nrow=nsim,ncol=2L-ncomp)
    routineName<-ifelse(ncomp==-1L,"mcMAM3","mcMAM4")
    ###
    res<-.Fortran(routineName,as.integer(nED),as.integer(nsim),as.double(ed1), 
                  as.double(sed1),as.double(addsigma),as.double(inis), 
                  as.integer(iflog),as.integer(nstart),as.double(w), 
                  as.double(m),chains=as.double(chains), 
                  iflag=as.integer(iflag),PACKAGE="numOSL")
    if (res$iflag!=0) {
        niter<-sum(res$chains[seq(nsim)]>0)
        stop(paste("Error: the simulation failed at the ", 
                   niter+1L,"th iteration!",sep=""))
    } # end if
    ###
    chains<-matrix(res$chains,ncol=2L-ncomp,
                   dimnames=list(NULL,
                   if(ncomp==-1L) c("p","gamma","sigma") else 
                   c("p","gamma","mu","sigma")))
    ###
    if (iflog==TRUE) {
        if (ncomp==-1L) {
            chains[,2L]<-exp(chains[,2L])
        } else if (ncomp==-2L) {
            chains[,2L:3L]<-exp(chains[,2L:3L])
        } # end if
    } # end if
    ###
    output<-list("EDdata"=EDdata, 
                 "addsigma"=addsigma, 
                 "model"=ifelse(ncomp==-1L,"MAM3","MAM4"),
                 "iflog"=iflog, 
                 "nsim"=nsim, 
                 "chains"=chains)
    class(output)<-"mcAgeModels"
    invisible(output)
    ###
} # end function mcMAM.
#####
