pattPC.fit<-function(obj, nitems,formel=~1,elim=~1,resptype="paircomp",obj.names=NULL,
         undec=TRUE, ia=FALSE, NItest=FALSE,NI=FALSE, MIScommon=FALSE,
         MISmodel="obj", MISalpha=NULL,MISbeta=NULL,
         pr.it=FALSE)
{
    call<-match.call()
    ENV<-new.env()
    ENV$pr.it<-pr.it

    ENV$resptype<-"paircomp"

    nobj<-nitems
    opt<-options()
    options("warn"=-1)
    ncomp<-nobj*(nobj-1)/2
    ENV$ncomp<-ncomp
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



   # for ratings: at least two items not NA (>1)
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

#######################################

    # check which response format (binary/ternary)
    ncat<-length(table(dat))
    if (ncat %% 2 > 0) {
         ncat<-3
    } else {
         ncat<-2
         undec<-FALSE
    }
    ENV$ncat <- ncat     # 7.12.09

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

    # missing values specifications
    #  check different combinations of MISalpha & MISbeta

    if(is.null(MISmodel)) MISmodel=" "

    if(is.null(MISalpha) || sum(MISalpha)==0 ){     # alpha not specified or all F
       if(is.null(MISbeta) || sum(MISbeta)==0 ){    #   beta not specified or all F
           MISalpha<-FALSE
           MISbeta<-FALSE
       } else {                                     #   beta specified
           ENV$Malph<-ENV$Mbeta<-MISbeta            #     alpha same spec as beta
           MISalpha<-TRUE
           MISbeta<-TRUE
       }
    } else {                                        # alpha specified
       if(is.null(MISbeta) || sum(MISbeta)==0 ){    #   beta not specified or all F
           MISbeta<-FALSE
           ENV$Malph<-MISalpha
           MISalpha<-TRUE
       } else {                                     #   beta specified
           ENV$Malph<-MISalpha
           ENV$Mbeta<-MISbeta
           MISalpha<-TRUE
           MISbeta<-TRUE
       }
    }
    # NI now parameter in function call
    if (NItest || MISalpha || MISbeta || MIScommon) NI <- TRUE # any nonresponse model treatment if TRUE
    ENV$NI <- NI
    if (MISalpha || MISbeta || MIScommon) NItest <- FALSE # NItest not useful with missing models
    if(MISbeta) MISalpha<-TRUE               # betas without alphas make no sense
    if(MIScommon) MISalpha<-MISbeta<-FALSE   # betas for common alpha make no sense

    ENV$NItest <- NItest
    ENV$MISalpha <- MISalpha
    ENV$MISbeta <- MISbeta
    ENV$MIScommon <- MIScommon

    if(is.null(MISmodel)) MISmodel=" "
    if(MISalpha && !(MISmodel %in% c("obj","comp")))
       #if(MISmodel %in% c("obj","comp"))
       #   ENV$MISmod <- MISmodel
       #else
          stop('\nMISmodel not correctly specified. Use "obj" or "comp"\n')
    ENV$MISmod <- MISmodel


    if(NI)
      if(!any(is.na(dat)))
         stop("no NA values - no nonresponse models !")

    if (ncat==2) {
          ENV$Y <- -(patternmat2(nobj)-1)            # pattern matrix
          ENV$Y<-ifelse(ENV$Y==1,1,-1)
    } else {
          ENV$Y <- -(patternmat3(nobj)-1)            # pattern matrix
    }

    if(undec)
      ENV$U <- apply(ENV$Y,1,function(x) sum(x==0))
    ENV$undec<-undec


    if(ENV$NI) {
          if(formel!="~1" || elim != "~1"){
          covs<-NULL
          formel<-~1
          elim<-~1
          cat("\ncurrently covariates fitted only for MCAR models !!\n")
       }
    }

    ENV$ia<-ia                            # dependency design matrix
    if (ia) {
       depL<-dependencies(nobj,ENV$Y)
       ENV$XI<-depL$d
       ilabels<-depL$label.intpars
       npars.ia<-nobj*(nobj-1)*(nobj-2)/2
       #npars.ia<-length(ilabels)
    } else {
       ilabels<-NULL
       npars.ia<-0
    }

    X<- ENV$Y %*% pcdesign(nobj)
    X <- X[,-nobj]                          # basic design matrix

    # split data according to subject covariates
    cList<-splitCovs(dat,covs,formel,elim,ENV)

    # generate list for all subj covariate x miss values groups
    if (ncat==2) {
          partsList<-gen.partsList2(nobj,cList,ENV)
    } else {
          partsList<-gen.partsList3(nobj,cList,ENV)
    }
    rm(cList)

    # number of parameters
        npar <- (nobj-1) * ENV$ncovpar + ENV$undec + npars.ia + NItest*(nobj-1)

    ########################################################
    #MNAR model Brian B:
    ##npar<-npar + MISalpha*(nobj) + MISbeta*(nobj)
    npar<-npar + MISalpha*sum(ENV$Malph) + MISbeta*sum(ENV$Mbeta) + MIScommon*1
    # 7.12.09 undec added to paridx
    paridx<-c(rep(1,nobj-1),rep(2,nobj*MISalpha)[ENV$Malph],rep(3,nobj*MISbeta)[ENV$Mbeta],4*MIScommon,6*undec,rep(5,ia*npars.ia))
    paridx<-paridx[paridx>0]
    ENV$paridx<-paridx
    ########################################################

    lambda<-rep(0.0,npar)
    ENV$iter<-0

    ## MAIN FITTING ROUTINE
    result<-nlm(loglikPC,lambda,X,nobj,partsList,ENV,hessian=TRUE,
        iterlim=1000)

#    prPC(result,ENV$covdesmat,nobj,elim,ilabels,ENV)
#    flush.console()
    options(opt)

    ENV$nobj<-nobj
    ENV$ilabels<-ilabels

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
