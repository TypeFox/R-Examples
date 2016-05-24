######################################################################
#
# rem.dyad.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 3/08/15
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# Contents:
#
#  accum.interact
#  acl.adj
#  acl.adjmat
#  acl.deg
#  accum.ps
#  acl.ps
#  acl.tri
#  accum.rrl
#  rem.dyad.lprior
#  rem.dyad.nlpost
#  rem.dyad.n2llik
#  rem.dyad.n2llik.samp
#  rem.dyad.nlpost.samp
#  rem.dyad.gof
#  rem.dyad
#  print.rem.dyad
#  print.summary.rem.dyad
#  summary.rem.dyad
#
######################################################################


#elist must be in the form cbind(time, src, dest), and time should be in
#ascending order.  Return value is a list with one element per event, 
#giving the accumulated activity level at that point.  Said levels are
#stored in lists of vectors, having length n.
accum.interact<-function(elist){
  acl<-.Call("accum_interact_R",as.matrix(elist),PACKAGE="relevent")
  names(acl)<-elist[,1]
  acl
}


#Return accumulated adjacency from src to dest by iteration iter
acl.adj<-function(acl,iter,src,dest){
  if(!(src%in%names(acl[[iter]])))
    return(0)
  else if(!(dest%in%names(acl[[iter]][[as.character(src)]])))
    return(0)
  else
    return(acl[[iter]][[as.character(src)]][[as.character(dest)]])
}
acl.adjmat<-function(acl,n,iter){
  a<-matrix(0,n,n)
  if(length(acl[[iter]])>0){
    for(i in 1:length(acl[[iter]])){
      for(j in 1:length(acl[[iter]][[i]])){
        a[as.numeric(names(acl[[iter]])[i]),as.numeric(names(acl[[iter]][[i]])[j])]<- acl[[iter]][[i]][[j]]
      }        
    }
  }
  a
}


#Return matrix of accumulated degrees
acl.deg<-function(acl,n,cmode=c("in","out","total")){
  if(match.arg(cmode)=="in"){
    deg<-matrix(0,length(acl),n)
    for(i in 1:length(acl))
      if(length(acl[[i]])>0){
        for(j in 1:length(acl[[i]])){  #Walk through persons with out-edges
          temp<-unlist(acl[[i]][[j]])  #Get tail names/values
          deg[i,as.numeric(names(temp))]<-deg[i,as.numeric(names(temp))]+temp
        }
      }
  }else if(match.arg(cmode)=="out"){
    deg<-matrix(0,length(acl),n)
    for(i in 1:length(acl))
      if(length(acl[[i]])>0){
        temp<-sapply(lapply(acl[[i]],unlist),sum) #Get outdegrees
        deg[i,as.numeric(names(acl[[i]]))]<-temp
      }
  }else{
    deg<-acl.deg(acl=acl,n=n,cmode="in")+acl.deg(acl=acl,n=n,cmode="out")
  }
  deg
}


#Build list of accumulated participation shift (P-shift) counts.  Output is in
#the form of a matrix whose rows contain counts for each of the 13 P-shift 
#types, such that the ith row contains the counts immediately prior to the
#resolution of event i.  Thus, there are m+1 rows in all, where m is the number
#of events.  The P-shifts are given in the order used in Gibson's 2003 Social
#Forces paper, namely:
#
#  (Turn Receiving)  [0] AB->BA, [1] AB->B0, [2] AB->BY,
#  (Turn Claiming)   [3] A0->X0, [4] A0->XA, [5] A0->XY,
#  (Turn Usurping)   [6] AB->X0, [7] AB->XA, [8] AB->XB, [9] AB->XY,
#  (Turn Continuing) [10] A0->AY, [11] AB->A0, [12] AB->AY  
#
#(This uses Gibson's notation, in which A is the initial source, B is the
#initial target, X is a new (shifted) speaker, Y is a new (shifted) target,
#and 0 is used where no well-defined speaker or target is present.  Here, this
#occurs when NA is given for source or destination.)
#
#It is worth noting that not all adjacent event pairs induce P-shifts, and hence
#the shift counts will not increment with every event.  In particular, the first
#event does not induce a shift (since there is no prior event), and neither does
#a repetition of a previous event (e.g., AB->AB or A0->A0).  The full set is
#thus affinely independent in general, although they will have a near (or 
#even full) dimension of affine dependence on most data sets.  You probably
#don't want to use all the stats at once, although we compute them all (since
#the cost of doing so is trivial).
accum.ps<-function(elist){
  psmat<-matrix(.Call("accum_ps_R",elist,PACKAGE="relevent"),ncol=13)
  colnames(psmat)<-c("AB-BA","AB-B0","AB-BY","A0-X0","A0-XA","A0-XY","AB-X0", "AB-XA","AB-XB","AB-XY","A0-AY","AB-A0","AB-AY")
  rownames(psmat)<-0:NROW(elist)
  psmat
}


#Return a list of 13 acl-style lists, each of which contains the changescores
#for one of the P-shift types.  The list structure is:
#  shifttype
#    $iter
#      $ego
#        $alter
#          $count (always 1)
acl.ps<-function(elist){
  n<-max(elist[,2:3],na.rm=TRUE)
  .Call("acl_ps_R",elist,n,PACKAGE="relevent")
}


#Return a list of pseudo-adjacencies reflecting triadic (really, 2-path) 
#properties:
# $ "top"|"tip"|"sop"|"sip"
#                          $iter
#                                $ ego
#                                      $ alter
#                                              $ count
acl.tri<-function(acl){
  .Call("acl_tri_R",acl,PACKAGE="relevent")
}


#Build list of accumulated incoming and outgoing recency-ranked communications.
#Specifically, the output is a nested list whose first dimension is in vs out,
#second dimension is iteration, third dimension is vertex, and fourth
#dimension (where applicable) is a recency ordered vector of communication
#partners (most recent to least recent).  This structure is built from
#the edgelist matrix, rather than the acl, due to the fact that this task
#is much easier to perform for the former than the latter.
accum.rrl<-function(elist){
  .Call("accum_rrl_R",as.matrix(elist),PACKAGE="relevent")
}


#Log prior density for the dyadic relational event model (assumes 
#multivariate t distribution), but will use a Gaussian limit at nu=Inf)
rem.dyad.lprior<-function(pv,pr.mean,pr.scale,pr.nu,...){
  if(all(pr.nu==Inf))
    sum(dnorm((pv-pr.mean)/pr.scale,log=TRUE)-log(pr.scale))
  else
    sum(dt((pv-pr.mean)/pr.scale,df=pr.nu,log=TRUE)-log(pr.scale))
}


#Negative log posterior for the dyadic relational event model
rem.dyad.nlpost<-function(pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar,ps,tri,lrm,pr.mean,pr.scale,pr.nu,ordinal,condnum,...){
  .Call("drem_n2llik_R",pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar, ps,tri,lrm,ordinal,condnum,PACKAGE="relevent")/2-rem.dyad.lprior(pv,pr.mean,pr.scale,pr.nu,...)
}


#Negative deviance for the dyadic relational event model
rem.dyad.n2llik<-function(pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar,ps,tri,lrm,ordinal,condnum,...){
  .Call("drem_n2llik_R",pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar, ps,tri,lrm,ordinal,condnum,PACKAGE="relevent")
}


#Negative deviance for the dyadic relational event model, using dyad sampling
rem.dyad.n2llik.samp<-function(pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar,ps,tri,lrv,tail,head,ordinal,condnum,...){
  .Call("drem_n2llik_samp_R",pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar, ps,tri,lrv,tail,head,ordinal,condnum,PACKAGE="relevent")
}


#Negative log posterior for the dyadic relational event model, using dyad 
#sampling
rem.dyad.nlpost.samp<-function(pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar,ps,tri,lrv,tail,head,pr.mean,pr.scale,pr.nu,ordinal,condnum,...){
  .Call("drem_n2llik_samp_R",pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl, covar,ps,tri,lrv,tail,head,ordinal,condnum,PACKAGE="relevent")/2-rem.dyad.lprior(pv,pr.mean,pr.scale,pr.nu,...)
}


#Goodness of fit statistics for the dyadic relational event model
rem.dyad.gof<-function(pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar,ps,tri,lrm,ordinal,condnum){
  .Call("drem_gof_R",pv,effects,edgelist,n,acl,cumideg,cumodeg,rrl,covar, ps,tri,lrm,ordinal,condnum,PACKAGE="relevent")
}


#Fit the dyadic relational event model, using the specified method and effects
rem.dyad<-function(edgelist,n,effects=NULL,ordinal=TRUE,acl=NULL,cumideg=NULL,cumodeg=NULL,rrl=NULL,covar=NULL,ps=NULL,tri=NULL,optim.method="BFGS",optim.control=list(),coef.seed=NULL,hessian=FALSE,sample.size=Inf,verbose=TRUE,fit.method=c("BPM","MLE","BSIR"),conditioned.obs=0,prior.mean=0,prior.scale=100,prior.nu=4,sir.draws=500,sir.expand=10,sir.nu=4,gof=TRUE){
  #t density
  dlmvt<-function(x,mu,is,ds,df){
    m<-length(x)
    lgamma((df+m)/2)-m/2*log(pi*df)-lgamma(df/2)-log(abs(ds))/2- (df+m)/2*log(1+((x-mu)%*%is%*%(x-mu))/df)
  }
#define NIDEGSEND  0   /*Norm indegree -> future sending rate*/
#define NIDEGREC   1   /*Norm indegree -> future receiving rate*/
#define NODEGSEND  2   /*Norm outdegree -> future sending rate*/
#define NODEGREC   3   /*Norm outdegree -> future receiving rate*/
#define NTDEGSEND  4   /*Norm total degree -> future sending rate*/
#define NTDEGREC   5   /*Norm total degree -> future receiving rate*/
#define FPSENDSEND 6   /*Fraction past sending -> future sending rate*/
#define FPRECSEND  7   /*Fraction past receipt -> future sending rate*/
#define RRRECSEND  8   /*Recency of receipt -> future sending rate*/
#define RRSENDSEND 9   /*Recency of sending -> future sending rate*/
#define COVSEND    10  /*Covariate effect for sending*/
#define COVREC     11  /*Covariate effect for receiving*/
#define COVSENDREC 12  /*Covariate effect for sending and receiving*/
#define COVEVENT   13  /*Generic event-wise covariate effect*/
#define OTPSEND    14  /*Outbound two-paths -> future sending rate*/
#define ITPSEND    15  /*Incoming two-paths -> future sending rate*/
#define OSPSEND    16  /*Outbound shared partners -> future sending rate*/
#define ISPSEND    17  /*Inbound shared partners -> future sending rate*/
#define FESEND     18  /*Fixed effects for sending*/
#define FEREC      19  /*Fixed effects for receiving*/
#define FESENDREC  20  /*Fixed effects for sending and receiving*/
#define PSABBA     21  /*P-Shift (turn receiving): AB->BA (dyadic)*/
#define PSABB0     22  /*P-Shift (turn receiving): AB->B0 (non-dyadic)*/
#define PSABBY     23  /*P-Shift (turn receiving): AB->BY (dyadic)*/
#define PSA0X0     24  /*P-Shift (turn claiming): A0->X0 (non-dyadic)*/
#define PSA0XA     25  /*P-Shift (turn claiming): A0->XA (non-dyadic)*/
#define PSA0XY     26  /*P-Shift (turn claiming): A0->XY (non-dyadic)*/
#define PSABX0     27  /*P-Shift (turn usurping): AB->X0 (non-dyadic)*/
#define PSABXA     28  /*P-Shift (turn usurping): AB->XA (dyadic)*/
#define PSABXB     29  /*P-Shift (turn usurping): AB->XB (dyadic)*/
#define PSABXY     30  /*P-Shift (turn usurping): AB->XY (dyadic)*/
#define PSA0AY     31  /*P-Shift (turn continuing): A0->AY (non-dyadic)*/
#define PSABA0     32  /*P-Shift (turn continuing): AB->A0 (non-dyadic)*/
#define PSABAY     33  /*P-Shift (turn continuing): AB->AY (dyadic)*/
  #Get effects
  allenam<-c("NIDSnd","NIDRec","NODSnd","NODRec","NTDegSnd","NTDegRec", "FrPSndSnd","FrRecSnd","RRecSnd","RSndSnd","CovSnd","CovRec","CovInt","CovEvent","OTPSnd","ITPSnd","OSPSnd","ISPSnd","FESnd","FERec","FEInt","PSAB-BA","PSAB-B0","PSAB-BY","PSA0-X0","PSA0-XA","PSA0-XY","PSAB-X0","PSAB-XA","PSAB-XB","PSAB-XY","PSA0-AY","PSAB-A0","PSAB-AY")
  effects<-allenam%in%effects
  #Fix the edgelist
  if(!is.matrix(edgelist))
    edgelist<-as.matrix(edgelist)
  if(ordinal&&is.na(edgelist[NROW(edgelist),2])) #Check for timestamp row
    edgelist<-edgelist[1:(NROW(edgelist)-1),,drop=FALSE]
  if(any(is.na(edgelist[1:(NROW(edgelist)-1+ordinal),]))){
    warning("Edgelist contains missing data; dropping incomplete events.\n")
    sel<-apply(is.na(edgelist),1,any)
    if(!ordinal)                     #Can allow NAs for timestamp row
      sel[length(sel)]<-TRUE
    edgelist<-edgelist[sel,,drop=FALSE]
  }
  if(any(diff(edgelist[,1])<=0,na.rm=TRUE))
    stop("Events are not well-ordered.  Stopping.\n")
  if(any(edgelist[,2]==edgelist[,3],na.rm=TRUE)){
    warning("Edgelist list contains loops (not currently supported); dropping self-interactions.\n")
    sel<-edgelist[,2]!=edgelist[,3]   #Get rid of loops...
    sel[is.na(sel)]<-TRUE
    edgelist<-edgelist[sel,]
  }
  #Check for covariates
  ncov<-rep(0,4)
  ctimed<-rep(FALSE,4)
  if(effects[11]){            #Covariate effects for sender
    if(is.null(covar$CovSnd))
      stop("CovSnd requires a corresponding covariate.\n")
    else{
      if(length(dim(covar$CovSnd))%in%(0:1)){
        ncov[1]<-1
        covar$CovSnd<-array(covar$CovSnd,dim=c(1,1,n))
        if(storage.mode(covar$CovSnd)!="double")
          storage.mode(covar$CovSnd)<-"double"
      }else if(length(dim(covar$CovSnd))==2){
        ncov[1]<-NCOL(covar$CovSnd)
        temp<-covar$CovSnd
        covar$CovSnd<-array(dim=c(1,ncov[1],n))
        covar$CovSnd[1,,]<-t(temp)
        if(storage.mode(covar$CovSnd)!="double")
          storage.mode(covar$CovSnd)<-"double"
      }else if(length(dim(covar$CovSnd))==3){
        ctimed[1]<-TRUE
        ncov[1]<-dim(covar$CovSnd)[2]
        if(dim(covar$CovSnd)[1]!=NROW(edgelist))
          stop("CovSnd covariate has ",dim(covar$CovSnd)[1]," time points, but needs",NROW(edgelist),".\n",sep="")
        if(storage.mode(covar$CovSnd)!="double")
          storage.mode(covar$CovSnd)<-"double"
      }else
        stop("CovSnd covariate has too many dimensions.\n")
    }
  }
  if(effects[12]){            #Covariate effects for receiver
    if(is.null(covar$CovRec))
      stop("CovRec requires a corresponding covariate.\n")
    else{
      if(length(dim(covar$CovRec))%in%(0:1)){
        ncov[2]<-1
        covar$CovRec<-array(covar$CovRec,dim=c(1,1,n))
        storage.mode(covar$CovRec)<-"double"
      }else if(length(dim(covar$CovRec))==2){
        ncov[2]<-NCOL(covar$CovRec)
        temp<-covar$CovRec
        covar$CovRec<-array(dim=c(1,ncov[2],n))
        covar$CovRec[1,,]<-t(temp)
        if(storage.mode(covar$CovRec)!="double")
         storage.mode(covar$CovRec)<-"double"
      }else if(length(dim(covar$CovRec))==3){
        ctimed[2]<-TRUE
        ncov[2]<-dim(covar$CovRec)[2]
        if(dim(covar$CovRec)[1]!=NROW(edgelist))
          stop("CovRec covariate has ",dim(covar$CovRec)[1]," time points, but needs",NROW(edgelist),".\n",sep="")
        if(storage.mode(covar$CovRec)!="double")
          storage.mode(covar$CovRec)<-"double"
      }else
        stop("CovRec covariate has too many dimensions.\n")
    }
  }
  if(effects[13]){            #Covariate effects for interaction
    if(is.null(covar$CovInt))
      stop("CovInt requires a corresponding covariate.\n")
    else{
      if(length(dim(covar$CovInt))%in%(0:1)){
        ncov[3]<-1
        covar$CovInt<-array(covar$CovInt,dim=c(1,1,n))
        if(storage.mode(covar$CovInt)!="double")
          storage.mode(covar$CovInt)<-"double"
      }else if(length(dim(covar$CovInt))==2){
        ncov[3]<-NCOL(covar$CovInt)
        temp<-covar$CovInt
        covar$CovInt<-array(dim=c(1,ncov[3],n))
        covar$CovInt[1,,]<-t(temp)
        if(storage.mode(covar$CovInt)!="double")
          storage.mode(covar$CovInt)<-"double"
      }else if(length(dim(covar$CovInt))==3){
        ctimed[3]<-TRUE
        ncov[3]<-dim(covar$CovInt)[2]
        if(dim(covar$CovInt)[1]!=NROW(edgelist))
          stop("CovInt covariate has ",dim(covar$CovInt)[1]," time points, but needs",NROW(edgelist),".\n",sep="")
        if(storage.mode(covar$CovInt)!="double")
          storage.mode(covar$CovInt)<-"double"
      }else
        stop("CovInt covariate has too many dimensions.\n")
    }
  }
  if(effects[14]){          #Covariate effects for events
    if(is.null(covar$CovEvent))
      stop("CovEvent requires a corresponding covariate.\n")
    else{
      if(length(dim(covar$CovEvent))%in%(0:1)){
        stop("CovEvent covariate has too few dimensions.\n")
      }else if(length(dim(covar$CovEvent))==2){
        ncov[4]<-1
        temp<-covar$CovEvent
        covar$CovEvent<-array(dim=c(1,ncov[4],n,n))
        covar$CovEvent[1,1,,]<-temp
        if(storage.mode(covar$CovEvent)!="double")
          storage.mode(covar$CovEvent)<-"double"
      }else if(length(dim(covar$CovEvent))==3){
        ncov[4]<-dim(covar$CovEvent)[1]
        covar$CovEvent<-array(covar$CovEvent,dim=c(1,ncov[4],n,n))
        if(storage.mode(covar$CovEvent)!="double")
          storage.mode(covar$CovEvent)<-"double"
      }else if(length(dim(covar$CovEvent))==4){
        ctimed[4]<-TRUE
        ncov[4]<-dim(covar$CovEvent)[2]
        if(dim(covar$CovEvent)[1]!=NROW(edgelist))
          stop("CovEvent covariate has ",dim(covar$CovEvent)[1]," time points, but needs",NROW(edgelist),".\n",sep="")
        if(storage.mode(covar$CovEvent)!="double")
          storage.mode(covar$CovEvent)<-"double"
      }else
        stop("CovEvent covariate has too many dimensions.\n")
    }
  }
  if(any(ncov>0)){
    covar$ncov<-as.integer(ncov)
    covar$ctimed<-as.logical(ctimed)
  }
  #Compute null deviance
  if(ordinal){
    nulldev<--2*(NROW(edgelist)-conditioned.obs)*log(1/(n*(n-1)))
    df.null<-0
  }else{
    timeint<-max(edgelist[,1])
    ecount<-NROW(edgelist)-1
    ecnt<-n*(n-1)
    erate<-ecount/timeint
    if(conditioned.obs==0){
      dt<-timeint
    }else{
      dt<-edgelist[ecount+1,1]-edgelist[conditioned.obs,1]
    }
    nulldev<- -2*(log(erate/ecnt)*(ecount-conditioned.obs)-erate*dt)
    df.null<-1
  }
  #Fit models, as needed
  if(sum(effects)==0){   #Return the null model
    fit<-list(n=n,m=NROW(edgelist)-conditioned.obs,ordinal=ordinal,null.deviance=nulldev,residual.deviance=nulldev,df.model=df.null, model.deviance=0,df.null=df.null,effects=effects,AIC=nulldev+2*df.null,AICC=nulldev+2*df.null*(df.null+1)/(NROW(edgelist)-conditioned.obs-(!ordinal)-df.null-1),BIC=nulldev+log(NROW(edgelist)-conditioned.obs-(!ordinal))*df.null)
  }else{
    #Perform initial setup
    stime<-proc.time()[3]
    if(verbose)
      cat("Computing preliminary statistics\n")
    if(is.null(acl)&&any(effects[1:8]))
      acl<-accum.interact(edgelist)
    if(is.null(cumideg)&&any(effects[1:8]))
      cumideg<-acl.deg(acl,n,cmode="in")
    if(is.null(cumodeg)&&any(effects[1:8]))
      cumodeg<-acl.deg(acl,n,cmode="out")
    if(is.null(rrl)&&any(effects[9:10]))
      rrl<-accum.rrl(edgelist)
    if(is.null(ps)&&any(effects[22:34]))
      ps<-acl.ps(edgelist)
    if(is.null(tri)&&any(effects[15:18])){
      if(is.null(acl))
        acl<-accum.interact(edgelist)
      tri<-acl.tri(acl)
    }
    if(sample.size<n*(n-1)){
      lrv<-rep(0.0,sample.size+1)
      temp<-rgnm(1,n,sample.size)
      tail<-as.integer(c(0,row(temp)[temp>0]-1))
      head<-as.integer(c(0,col(temp)[temp>0]-1))
      lrm<-NULL
    }else{
      lrm<-matrix(0.0,n,n)
      lrv<-NULL
      head<-NULL
      tail<-NULL
    }
    nparm<-sum(effects[c(1:10,15:18,22:34)])+sum(ncov)+sum(effects[19:21]*(n-1))
    if(is.null(coef.seed))
      pv<-rnorm(nparm,0,0.001)
    else
      pv<-coef.seed
    if(match.arg(fit.method)!="MLE"){
      prior.mean<-rep(prior.mean,length=nparm)
      prior.scale<-rep(prior.scale,length=nparm)
      prior.nu<-rep(prior.nu,length=nparm)
    }
    #Fit the model
    if(verbose)
      cat("Fitting model\n")
    if(sample.size<n*(n-1)){
      if(match.arg(fit.method)=="MLE"){
        ofun<-match.fun("rem.dyad.n2llik.samp")
      }else{
        ofun<-match.fun("rem.dyad.nlpost.samp")
        lfun<-match.fun("rem.dyad.n2llik.samp")
      }
    }else{
      if(match.arg(fit.method)=="MLE"){
        ofun<-match.fun("rem.dyad.n2llik")
      }else{
        ofun<-match.fun("rem.dyad.nlpost")
        lfun<-match.fun("rem.dyad.n2llik")
      }
    }
    if(match.arg(fit.method)=="BSIR")
      hessian<-TRUE
    fit<-optim(pv,ofun,method=optim.method,control=optim.control, hessian=hessian,effects=effects,edgelist=edgelist,n=n,acl=acl,cumideg=cumideg, cumodeg=cumodeg,rrl=rrl,covar=covar,ps=ps,tri=tri,lrv=lrv,lrm=lrm,tail=tail,head=head,pr.mean=prior.mean,pr.scale=prior.scale,pr.nu=prior.nu,ordinal=ordinal,condnum=conditioned.obs)
    fit$coef<-fit$par
    cnam<-c(allenam[1:10],paste(allenam[11],1:(max(1,ncov[1])),sep="."), paste(allenam[12],1:(max(1,ncov[2])),sep="."), paste(allenam[13],1:(max(1,ncov[3])),sep="."), paste(allenam[14],1:(max(1,ncov[4])),sep="."),allenam[15:18], paste(allenam[19],2:n,sep="."),paste(allenam[20],2:n,sep="."), paste(allenam[21],2:n,sep="."),allenam[22:34])[c(effects[1:10], rep(effects[11:14],times=pmax(1,ncov)),effects[15:18], rep(effects[19:21],each=n-1),effects[22:34])]
    names(fit$coef)<-cnam
    fit$n<-n
    fit$m<-NROW(edgelist)-conditioned.obs
    fit$ordinal<-ordinal
    fit$conditioned.obs<-conditioned.obs
    if(gof){
      if(verbose)
        cat("Obtaining goodness-of-fit statistics\n")
      if(sample.size<n*(n-1))
        lrm<-matrix(0.0,n,n)
      temp<-rem.dyad.gof(fit$par,effects=effects,edgelist=edgelist, n=n,acl=acl,cumideg=cumideg,cumodeg=cumodeg,rrl=rrl,covar=covar,ps=ps,tri=tri,lrm=lrm,ordinal=ordinal,condnum=conditioned.obs)
      fit$residual.deviance<-sum(temp$residuals)+temp$dev.censor
      fit$model.deviance<-nulldev-fit$residual.deviance
      fit$residuals<-temp$residuals
      fit$predicted<-matrix(temp$predicted,ncol=2)
      fit$predicted.match<-temp$predicted==edgelist[(1+conditioned.obs): (NROW(edgelist)-1+ordinal), 2:3]
      fit$observed.rank<-temp$obs.rank
    }else{
      fit$residual.deviance<-fit$value #Deviance if MLE, else -log post
      if(match.arg(fit.method)%in%c("BPM","BSIR")){  #Correct to dev if not MLE
        fit$residual.deviance<-2*(fit$residual.deviance + rem.dyad.lprior(pv=fit$coef,pr.mean=prior.mean,pr.scale=prior.scale,pr.nu=prior.nu))
      }
      fit$model.deviance<-nulldev-fit$residual.deviance
    }
    fit$df.model<-nparm
    fit$null.deviance<-nulldev
    fit$df.null<-df.null
    fit$effects<-effects
    fit$AIC<-fit$residual.deviance+2*nparm
    fit$AICC<-fit$AIC+2*nparm*(nparm+1)/(fit$m-(!ordinal)-nparm-1)
    fit$BIC<-fit$residual.deviance+nparm*log(fit$m-(!ordinal))
    if(!is.null(fit$hessian))
      fit$cov<-qr.solve(fit$hessian)
    if(match.arg(fit.method)=="BSIR"){
      if(verbose)
        cat("Taking posterior draws\n")
      cf<-chol(fit$cov)
      is<-qr.solve(fit$cov)
      ds<-det(fit$cov)
      post<-matrix(0,sir.draws*sir.expand,nparm)
      iw<-vector()
      dev<-vector()
      lp<-vector()
      for(i in 1:(sir.draws*sir.expand)){
        post[i,]<-fit$coef+cf%*%rnorm(nparm)*sqrt(sir.nu/rchisq(1,sir.nu))
        dev[i]<-lfun(post[i,],effects=effects,edgelist=edgelist,n=n,acl=acl, cumideg=cumideg,cumodeg=cumodeg,rrl=rrl,covar=covar,ps=ps,tri=tri,lrv=lrv,lrm=lrm,tail=tail,head=head,ordinal=ordinal,condnum=conditioned.obs)
        lp[i]<-rem.dyad.lprior(post[i,],pr.mean=prior.mean,pr.scale=prior.scale, pr.nu=prior.nu)
        iw[i]<--dev[i]/2+lp[i]-dlmvt(post[i,],fit$coef,is,ds,sir.nu)
        #cat("post:\n\t")
        #print(post[i,])
        #cat("dev=",dev[i],"lp=",lp[i],"spr=",dlmvt(post[i,],fit$coef,is,ds,3), "iw=",iw[i],"\n")
      }
      iw<-iw-logSum(iw)
      if(verbose){
        if(any(is.na(iw)|is.nan(iw)|is.na(exp(iw)))||all(exp(iw)==0)){
          print(cbind(post,dev,lp,iw,exp(iw)))
        }
        print(quantile(iw,(0:10)/10))
        print(quantile(exp(iw),(0:10)/10))
        #hist(exp(iw))
      }
      sel<-sample(1:NROW(post),sir.draws,replace=TRUE,prob=exp(iw))
      fit$post<-post[sel,]
      fit$iw<-iw
      fit$post.deviance<-dev[sel]
      fit$post.mean.deviance<-mean(dev[sel])
      fit$DIC<-2*fit$post.mean.deviance-fit$residual.deviance
      fit$pd1<-fit$post.mean.deviance-fit$residual.deviance
      fit$pd2<-sd(dev[sel])
      fit$post.cov<-var(fit$post)
    }
    fit$compute.time<-proc.time()[3]-stime
  }
  class(fit)<-"rem.dyad"
  #Return the result
  fit
}

print.rem.dyad<-function(x, ...){
  cat("Relational Event Model\n")
  if(is.null(x$coef)){
    cat("\nNull model object.\n")
  }else{
    cat("\nCoefficient Estimates:\n")
    print(x$coef)
    #print(x$effects)
  }
  cat("\nNull Deviance:",x$null.deviance,"\n")
  cat("Residual Deviance:",x$residual.deviance,"\n")
  cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n\n")
}

print.summary.rem.dyad<-function(x, ...){
  cat("Relational Event Model ")
  if(x$ordinal)
    cat("(Ordinal Likelihood)\n\n")
  else
    cat("(Temporal Likelihood)\n\n")
  if(is.null(x$coef)){
    cat("Null model object.\n\n")
  }else{
    if(is.null(x$cov)){
      ctab<-matrix(x$coef,ncol=1)
      rownames(ctab)<-names(x$coef)
      colnames(ctab)<-c("Estimate")
      printCoefmat(ctab)
    }else{
      ctab<-cbind(x$coef,diag(x$cov)^0.5)
      ctab<-cbind(ctab,ctab[,1]/ctab[,2])
      ctab<-cbind(ctab,2*(1-pnorm(abs(ctab[,3]))))
      rownames(ctab)<-names(x$coef)
      colnames(ctab)<-c("Estimate","Std.Err","Z value","Pr(>|z|)")
      printCoefmat(ctab,P.values=TRUE)
    }
  }
  cat("Null deviance:",x$null.deviance,"on",x$m-x$df.null,"degrees of freedom\n")
  cat("Residual deviance:",x$residual.deviance,"on",x$m-x$df.model,"degrees of freedom\n")
  cat("\tChi-square:",x$model.deviance,"on",x$df.model-x$df.null,"degrees of freedom, asymptotic p-value",1-pchisq(x$model.deviance,x$df.model-x$df.null),"\n")
  cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n")
}

summary.rem.dyad<-function(object, ...){
  class(object)<-c("summary.rem.dyad",class(object))
  object
}


