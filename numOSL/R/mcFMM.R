#####
mcFMM<-
function(EDdata, ncomp=1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    UseMethod("mcFMM")
} #
### 2014.09.20.
mcFMM.default<-
function(EDdata, ncomp=1, addsigma=0, iflog=TRUE,
         nsim=5e4, inis=list(), control.args=list()) {
    ### Stop if not.    
     stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,2L,drop=TRUE]>0),
              length(ncomp)==1L, ncomp %in% seq(4L),
              length(addsigma)==1L, is.numeric(addsigma),
              length(iflog)==1L, is.logical(iflog),
              length(nsim)==1L, is.numeric(nsim), nsim>=100L, nsim<=2e5,
              is.list(inis), is.list(control.args),
              (ncomp==1L && all(names(inis) %in% c("mu","sigma"))) ||
              (ncomp>1L && all(names(inis) %in% c(paste("p",seq(ncomp),sep=""),
               paste("mu",seq(ncomp),sep="")))),
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
    ### Default initials.
    args.inis<-list("mu"=mean(ed1), 
                    "sigma"=0.5,
                    "p1"=1.0, "p2"=1.0,
                    "p3"=1.0, "p4"=1.0,
                    "mu1"=mean(ed1), 
                    "mu2"=mean(ed1),
                    "mu3"=mean(ed1),
                    "mu4"=mean(ed1))
    args.inis[names(inis)]<-inis
    stopifnot(args.inis[["mu"]]>lowerGamma,
              args.inis[["mu"]]<upperGamma,
              args.inis[["sigma"]]>0,
              args.inis[["sigma"]]<ifelse(iflog==TRUE,5,var(ed1)),
              args.inis[["p1"]]>0, args.inis[["p2"]]>0,
              args.inis[["p3"]]>0, args.inis[["p4"]]>0,
              args.inis[["mu1"]]>lowerGamma, 
              args.inis[["mu1"]]<upperGamma,
              args.inis[["mu2"]]>lowerGamma, 
              args.inis[["mu2"]]<upperGamma,
              args.inis[["mu3"]]>lowerGamma, 
              args.inis[["mu3"]]<upperGamma,
              args.inis[["mu4"]]>lowerGamma, 
              args.inis[["mu4"]]<upperGamma)
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
    inis<-if(ncomp==1L) {
        c(args.inis[["mu"]],args.inis[["sigma"]])
    } else if (ncomp==2L) {
        matrix(c(args.inis[["p1"]],args.inis[["mu1"]],
                 args.inis[["p2"]],args.inis[["mu2"]]),
                 ncol=2L)
    } else if (ncomp==3L) {
        matrix(c(args.inis[["p1"]],args.inis[["mu1"]],
                 args.inis[["p2"]],args.inis[["mu2"]],
                 args.inis[["p3"]],args.inis[["mu3"]]),
                 ncol=3L)
    } else if (ncomp==4L) {
        matrix(c(args.inis[["p1"]],args.inis[["mu1"]],
                 args.inis[["p2"]],args.inis[["mu2"]],
                 args.inis[["p3"]],args.inis[["mu3"]],
                 args.inis[["p4"]],args.inis[["mu4"]]),
                 ncol=4L)
    } # end if
    iflag<-0
    chains<-matrix(0,nrow=nsim,ncol=2L*ncomp)
    ###
    routineName<-ifelse(ncomp==1L, "mcCAM", 
                 paste("mcFMM",ncomp,sep=""))
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
    chains<-matrix(res$chains,ncol=2L*ncomp)
    if (ncomp==1L) {
        if (iflog==TRUE) {
            chains[,1L]<-exp(chains[,1L])
        } else {
            chains[,2L]<-chains[,2L]/chains[,1L]
        } # end if
    } else {
        if (iflog==TRUE) {
            chains[,(ncomp+1L):(2L*ncomp)]<-
            exp(chains[,(ncomp+1L):(2L*ncomp)])
        } # end if
        ### Constraint with u1<u2<u3<u4.
        chains<-t(apply(array(t(as.matrix(chains)),dim=c(ncomp,2L,nsim)), 
                  MARGIN=3L, function(x)  x[order(x[,2L,drop=TRUE]),]))
    } # end if
    ###
    dimnames(chains)=list(NULL,
                     if(ncomp==1L) c("mu","sigma") else 
                     c(paste("p",seq(ncomp),sep=""),
                     paste("mu",seq(ncomp),sep="")))
    ### 
    output<-list("EDdata"=EDdata, 
                 "addsigma"=addsigma, 
                 "model"=ifelse(ncomp==1L,"CAM",paste("FMM",ncomp,sep="")), 
                 "iflog"=iflog, 
                 "nsim"=nsim, 
                 "chains"=chains)
    class(output)<-"mcAgeModels"
    invisible(output)
    ###
} # end function mcFMM.
#####          
