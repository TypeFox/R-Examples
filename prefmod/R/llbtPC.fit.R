####### function for LLBT - from pattPC.fit   #############################################
llbtPC.fit<-function(obj, nitems, formel=~1, elim=~1, resptype="paircomp", obj.names=NULL, undec=TRUE)
{
    call<-match.call()
    ENV<-new.env()

    ENV$resptype<-"paircomp"
    ENV$undec <- undec # added 2011-11-01

    nobj<-nitems
    ENV$nobj<-nobj
    opt<-options()
    options("warn"=-1)
    ncomp<-nobj*(nobj-1)/2
    ENV$ncomp<-ncomp

    ## process inout data
    if(is.character(obj)){                    ## datafilename supplied
         datafile    <-  obj
         if(file.access(datafile, mode=0) == 0){
            dat<-as.matrix(read.table(datafile,header=TRUE))  # datafile
         } else {
            stop("\ninput data file does not exist!\n")
         }
    } else if(is.data.frame(obj)){            ## data frame supplied
         dat<-as.matrix(obj)                                  # dataframe
         dat<-apply(dat,2,as.numeric)
    } else {
         stop("first argument must be either datafilename or dataframe")
    }
    varnames<-colnames(dat)
    if (ncol(dat)>ncomp) {
         covnames<-varnames[(ncomp+1):ncol(dat)]
         covs<-as.data.frame(dat[,(ncomp+1):ncol(dat)])
    } else {
         covs<-NULL
    }
    # at least one comparison not NA
    idx<-apply(dat[,1:ncomp],1,function(x) sum(!is.na(x))>=1)
    dat<-dat[idx,]
    dat<-as.data.frame(dat[,1:ncomp])

    if(!is.null(covs)){
         covs<-as.data.frame(covs[idx,])
         colnames(covs)<-covnames
         NAs<-which(!complete.cases(covs))                  # check for NA
         if (length(NAs)>0){
              cat("\tsubject covariates: NAs in lines",NAs," - removed from data\n")
              notNAs<-which(complete.cases(covs))
              dat<-dat[notNAs,]
              # covs<-covs[notNAs,] ## replaced 20-08-09
              covs<-covs[notNAs,,drop=FALSE]
         }
    }

    # check which response format (binary/ternary)
    ncat<-length(table(as.matrix(dat)))
    if (ncat %% 2 > 0) {
         ncat<-3
         pcdes.mult <- 1:-1
    } else {
         ncat<-2
         pcdes.mult <- c(1,-1)
         undec<-FALSE
    }

    # transform data to 0,1 or 0,1,2
    if (min(dat, na.rm=TRUE)<0) dat <- -dat  # if -1/1 used then larger value is the preferred
    dat<-as.data.frame(red.cat(dat,ncat)) # changes data structure to 0,1 or 0,1,2

    # names for objects/items
    if (!is.null(obj.names)) {
        if (is.character(obj.names) && length(obj.names)==nobj) {
          ENV$obj.names<-obj.names
        } else {
          cat("\nobj.names inccorrectly specified - o1, o2, .. is used instead\n")
          ENV$obj.names<-paste("o",1:nobj,sep="")
        }
    } else {
          ENV$obj.names<-paste("o",1:nobj,sep="")
    }


    ## basic design matrix

    # split data according to subject covariates
    cList<-splitCovs(dat,covs,formel,elim,ENV)

    # calculate y
#    cLtab <- lapply(cList, function(x) apply(as.matrix(as.data.frame(x[[1]])),2,table))
    cLtab <- lapply(cList, function(x) lapply(x[[1]],function(y) tabulate(na.omit(y+1),nbins=ncat) ))
    y <- unlist(cLtab)

    # calculate mu
    mu <- gl(ncomp , ncat, ncat*ncomp*length(cList))

    dfr <- data.frame(y, mu)

    # undecided
    if(undec){
        U <- rep(c(0,1,0),length(y)/3)
        dfr<-cbind(dfr,U)
    }

    # objects design
    pcdes <- rep(1,length(cList)) %x% pcdesign(nobj)
    objdesign <- pcdes %x% pcdes.mult
    colnames(objdesign) <- ENV$obj.names

    dfr <- cbind(dfr, objdesign)
    # covariates
    if (length(cList[[1]]$cov) > 1) {

       covsdesign <- as.data.frame(gfac2(ENV$elimcovlevels) %x% rep(1,ncomp*ncat))
       names(covsdesign)<-names(ENV$elimcovlevels)
       covsdesign <- covsdesign[names(ENV$model.covs)]

       #covsdesign <- data.frame(gfac2(ENV$covlevels) %x% rep(1,ncomp*ncat)) # 23.11.09 if based on model terms only
       #names(covsdesign) <- names(ENV$covlevels)                            # we get a wrong sequence of factor levels

       names(covsdesign) <- names(ENV$model.covs)
       covsdesign <- lapply(covsdesign, function(x) as.factor(x))
       dfr <- cbind(dfr, covsdesign)
    }

 #ev. gnm.Fit


    ## model formula

    # part for objects
    frm.objects <- paste("(",paste(ENV$obj.names,sep="",collapse="+"),")",sep="",collapse="")

    # part for covs
    if (length(cList[[1]]$cov) > 1) {
       # rh 22.4.10: use ENV$formel instead of formel because ENV$formel has sorted terms
       frm.covs <- paste("+",frm.objects, ":(",gsub("[~[:blank:]]","",ENV$formel)[2],")",sep="",collapse="")
    } else {
       frm.covs <- ""
    }

    # part for undecided
    if (undec) {
       frm.u  <- "+ U"
    } else {
       frm.u <- ""
    }

    formula <- as.formula(paste("y~",frm.objects,frm.covs,frm.u, sep="",collapse=""))

    npar.elim <- ncomp * length(cList)
    elim <- gl(npar.elim, ncat)
    dfr <- cbind(dfr, elim)


    rm(cList,cLtab)
    gc()

    result <- gnm(formula, eliminate=elim, family=poisson, data=dfr)
    result$envList <- as.list(ENV)
    class(result) <- c("llbtMod", class(result))
    result
}
