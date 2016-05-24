pattLrep.fit<-function(obj, nitems, tpoints=1, formel=~1,elim=~1,resptype="ratingT",
         obj.names=NULL, undec=TRUE, ia=FALSE, iaT=FALSE, NItest=FALSE, pr.it=FALSE)
{
    # nitems ... no of items at a given time point
    # tpoints ... no of timepoints

    if (tpoints<2)
        stop("no of timepoints incorrectly specified! if tpoints==1 use pattL.fit")

    call<-match.call()
    ENV<-new.env()
    ENV$pr.it<-pr.it

    ENV$resptype<-"ratingT"

    nobj<-nitems * tpoints
    opt<-options()
    options("warn"=-1)

#######################################
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
   if (ncol(dat)>nobj) {                                         # from pattR.fit 2011-08-31
        #covnames<-varnames[(nobj+1):ncol(dat)]
        #covs<-as.data.frame(dat[,(nobj+1):ncol(dat)])
        ## instead of the above rh 2011-05-13
        formel.names<-attr(terms(as.formula(formel)),"term.labels")
        formel.names<-unique(unlist(strsplit(formel.names,":"))) # 2011-08-31 remove interaction terms
        elim.names<-attr(terms(as.formula(elim)),"term.labels")
        elim.names<-unique(unlist(strsplit(elim.names,":")))     # 2011-08-31 remove interaction terms
        covnames<-unique(c(formel.names,elim.names))
        covs<-as.data.frame(dat[,covnames])
   } else {
        covs<-NULL
   }

   # for ratings: at least two items not NA
   idx<-apply(dat[,1:nobj],1,function(x) sum(!is.na(x))>1)
   dat<-dat[idx,]
   dat<-as.data.frame(dat[,1:nobj])

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

#######################################

    ## option for obj.names added
    if (is.null(obj.names))
      ENV$obj.names<-varnames[1:nobj]
    else
      ENV$obj.names<-obj.names[1:nobj]

    if(NItest)
      if(!any(is.na(dat)))
         stop("Test for ignorable missing cannot be performed - no NA values!")

    datrng<-range(dat,na.rm=TRUE)

    ## reduced patternmatrix for all timepoints
    # pattern vary fastest in left columns (first timepoint)

    Y <- Lpatternmat(datrng,nitems) # first pattern matrix reduced

    # "recursively" expand patternmatrix for timepoints
    np<-nrow(Y)                     # no patterns in Y
    npp<-np
    YL<-Y
    for (t in 1:(tpoints-1)){
      YL<-do.call("rbind", lapply(1:np, function(i) YL) ) # stacks YL np times
      YR<-expand.mat(Y,rep(npp,np))                       # repeats each line of Y np^t times
      YL<-cbind(YL,YR)
      npp<-npp*np
    }
    ENV$Y<-YL
    rm(Y,YL) # tidy up

    ## reduce datamatrix
    # for each timepoint
    #dat.t<-NULL
    dat.t<-as.data.frame(diffsred(dat[,1:nitems],nitems))
    for (t in 2:tpoints){
       from<-nitems*(t-1)+1
       to<-from+nitems-1
       dat.t<-cbind(dat.t,as.data.frame(diffsred(dat[,from:to],nitems)))
    }
    dat<-dat.t
    rm(dat.t)


    ncomp<-choose(nitems,2)
    # only global undecided
    if(undec){
      ENV$U <- apply(ENV$Y[,1:ncomp],1,function(x) sum(x==0))
      for (t in 2:tpoints){
         from<-ncomp*(t-1)+1
         to<-from+ncomp-1
         ENV$U <- cbind(ENV$U,apply(ENV$Y[,from:to],1,function(x) sum(x==0)))
      }
    }
    ENV$undec<-undec

    ENV$NItest<-NItest
    if(ENV$NItest) {
          if(formel!="~1" || elim != "~1"){
          covs<-NULL
          formel<-~1
          elim<-~1
          cat("\ncurrently no covariates fitted if NItest==TRUE !!\n")
       }
    }

    # dependence parameters within each time point
    ENV$ia<-ia
    ilabels<-NULL
    XI<-NULL
    if (ia) {
       XI<-NULL
       ilabels<-NULL
       for (t in 1:tpoints){
          from<-ncomp*(t-1)+1
          to<-from+ncomp-1
          depL<-dependencies(nitems,ENV$Y[,from:to])
          XI<-cbind(XI,depL$d)
          ilabels<-c(ilabels,depL$label.intpars)
       }
       #npars.ia<-tpoints * nitems*(nitems-1)*(nitems-2)/2
       npars.ia<-length(ilabels)
       ilabels<-paste(rep(paste("T",1:tpoints,":",sep=""),
                            rep(npars.ia/tpoints,tpoints)),ilabels,sep="")
    } else {
       ENV$ilabels<-NULL
       npars.ia<-0
    }
    ENV$XI<-XI
    rm(XI)
    ENV$ilabels<-ilabels

    # dependence parameters between timepoints (AR(1))
    ENV$iaT<-iaT
    if (iaT) {
       npars.iaT<-ncomp*(tpoints-1)
       ENV$XIT<-do.call("cbind", lapply(1:npars.iaT,function(i) ENV$Y[,i]*ENV$Y[,i+ncomp]))
       ENV$iTlabels<-paste(paste("Comp",1:ncomp,sep=""),
                       paste("IT",rep(1:(tpoints-1),rep(ncomp,tpoints-1)),rep(2:(tpoints),rep(ncomp,tpoints-1)),sep=""),
                    sep=":")
    } else {
       ENV$iTlabels<-NULL
       npars.iaT<-0
    }

    ncomp<-choose(nitems,2)
    X<- -(ENV$Y[,1:ncomp] %*% pcdesign(nitems))[,-nitems]
    for (t in 2:tpoints){
       from<-ncomp*(t-1)+1
       to<-from+ncomp-1
       X<-cbind(X,-(ENV$Y[,from:to] %*% pcdesign(nitems))[,-nitems]  )
    }

    #X<- ENV$Y %*% pcdesign(nobj)
    #X<- -X[,-nobj]                          # basic design matrix

    cList<-splitCovs(dat,covs,formel,elim,ENV)   # split data according to subject covariates
    partsList<-gen.partsList(nobj,cList,ENV)     # generate list for all subj covariate x miss values groups
    rm(cList)


    npar <- tpoints*(nitems-1) * ENV$ncovpar + ENV$undec*tpoints + npars.ia + npars.iaT
    if (ENV$NItest) npar<-tpoints*(nitems-1)*2 + ENV$undec*tpoints + npars.ia + npars.iaT

    lambda<-rep(0,npar)
    ENV$iter<-0


    nobj<-tpoints*(nitems-1)

    ## MAIN FITTING ROUTINE
    result<-nlm(loglik,lambda,X,nobj,partsList,ENV,hessian=TRUE,
        iterlim=1000)

    if (pr.it) cat("\n")
    options(opt)

    ENV$nobj<-nobj
    ENV$nitems<-nitems
    ENV$tpoints<-tpoints

    envList<-mget(ls(ENV),envir=ENV)
    outputobj<-list(coefficients=result$estimate,
                    ll=ENV$ll,
                    fl=ENV$fl,
                    call=call,
                    result=result,
                    envList=envList,
                    partsList=partsList)
    class(outputobj) <- c("pattMod")                         #class: pattern model
    outputobj
}
