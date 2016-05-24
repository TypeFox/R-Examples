
#_______________________________________________________________________________
# MAPLES library
# this file contains the following objects:
# Functions: ageprofile,epdata, splitter,  plotap, tabx, tabm, mkdate, listvar
#                         
#_______________________________________________________________________________




#_______________________________________________________________________________
# ageprofile   estimates age profiles
#_______________________________________________________________________________

ageprofile<-function(formula,epdata,
        tr.name="Transition",win,method="car",
        agelimits=c(0,100),outfile=FALSE,
        tails=c(FALSE,FALSE),subset=TRUE,
        weight,sig.eff=TRUE, sig.lev=0.05)
{

  #_______________________________________________________________________________
  # aggreg   aggregates micro data into a transition-specific matrix
  aggreg<-function(e,method,covk)
  {
    ncat<-nlevels(factor(covk))  # number of levels of covariate
    categ<- levels(factor(covk))   # levels of covariate
    event<-matrix(rep(0,100*ncat),ncol=100) # initialize event counter
    expos<-matrix(rep(0,100*ncat),ncol=100) # initialize exposure counter

    # cohort-age rates
    if (method=="car"){
      for (k in 1:ncat)  {
        ds<-subset(e,covk==categ[k])#Selects temporary subset
        N<-dim(ds)[1]
        ti<-pmax(trunc(ds$Astart),0)
        tf<-pmin(trunc(ds$Astop),100)
        st<-ifelse(ds$status==1 | ds$status==2,1,0)
        if (N>0) {
           for (j in 1:N) {    # selects j-th individual
              expos[k,][ti[j]:tf[j]] <-expos[k,][ti[j]:tf[j]]+ ds$weight[j]
              expos[k,][tf[j]] <-expos[k,][tf[j]] -(1-(ds$Astop[j]-
                            tf[j]))*ds$weight[j]
              expos[k,][ti[j]]<-expos[k,][ti[j]]-((ds$Astart[j]-ti[j])*ds$weight[j])
              event[k,][tf[j]]<-event[k,][tf[j]]+st[j]*ds$weight[j]
    }}}}
    # cohort-period rates
    if (method=="cpr"){
      for (k in 1:ncat)  {
        ds<-subset(e,covk==categ[k]) #Selects temporary subset
        N<-dim(ds)[1]
        ti<-pmax(trunc(trunc(ds$Wstart)-ds$birth)+1,0)
        tf<-pmin(trunc(trunc(ds$Wstop)-ds$birth)+1,100)
        st<-ifelse(ds$status==1 | ds$status==2,1,0)
      if (N>0) {
        for (j in 1:N) {  # selects j-th individual
          expos[k,][ti[j]:tf[j]] <-expos[k,][ti[j]:tf[j]]+ ds$weight[j]
          expos[k,][tf[j]]<-expos[k,][tf[j]]-
                          ((1-(ds$Wstop[j]-trunc(ds$Wstop[j])))*ds$weight[j])
          expos[k,][ti[j]]<-expos[k,][ti[j]]-
                           (ds$Wstart[j]-trunc(ds$Wstart[j]))*ds$weight[j]
          event[k,][tf[j]]<-event[k,][tf[j]]+st[j]*ds$weight[j]
    }}}} 
    # Cohort-age probabilities
    if (method=="cap"){
      for (k in 1:ncat)  {
        ds<-subset(e,covk==categ[k]) #Selects temporary subset
        N<-dim(ds)[1]
        ti<-pmax(trunc(ds$Astart),0)
        tf<-pmin(trunc(ds$Astop),100)
        st<-ifelse(ds$status==1 | ds$status==2,1,0)
        if (N>0) {
            for (j in 1:N) {  # selects j-th individual
              expos[k,][ti[j]:tf[j]] <-expos[k,][ti[j]:tf[j]]+ ds$weight[j]
          event[k,][tf[j]]<-event[k,][tf[j]]+st[j]*ds$weight[j]
    }}}} 
  # Recomposes data matrix
  datr<-data.frame(matrix(rep(0,400*ncat),ncol=4))
  names(datr)<-c("age","event","expos","variab")
  for (i in 1:100) {
    for (k in 1:ncat) {
       datr[ncat*(i-1)+k,1]<-i           #age
       datr[ncat*(i-1)+k,2]<-event[k,i]  #events
       datr[ncat*(i-1)+k,3]<-expos[k,i] #exposures
       datr[ncat*(i-1)+k,4]<-categ[k]  # variab
    }
  }
  return(datr)
  }

  #_____________________________________________________________________________
  # tailflat   Tail flattening procedure (Baseline smoothing at the hedges)
  tailflat<-function(age,profile,prop,tails)
  {
    smoo<-profile
    if (tails[1]){
      lim<- min(age[which(prop>=0.10)])   
      fdist<-which(age==lim)-1
      fw<-c(rep(0,fdist))
      r<-10/fdist
      cc<-exp(((fdist+1)*r)/2)
      fweight<-1/(1+cc*exp(-r*seq(1:fdist)))
      if (lim>min(age)) smoo[age<lim[1]]<-profile[age<lim]*fweight
    }     
    if (tails[2]){
      lim<-min(age[which(prop>=0.90)])
      fdist<-which(age==max(age))-which(age==lim)
      fw<-c(rep(0,fdist))
      r<-10/fdist
      cc<-exp(((fdist+1)*r)/2)
      fweight<-1/(1+cc*exp(-r*seq(1:fdist)))
      if (lim<max(age)) smoo[age>lim]<-profile[age>lim]*(1-fweight)
    }
    return(smoo)
  }  
  
  
  # CHECK Arguments ____________________________________________________________
  if (missing(formula)) stop ("Formula is missing (try with formula=~1)"
    ,call.=FALSE)
  if (missing(epdata)) stop("Must have an epdata argument",call.=FALSE)
  if (!is.data.frame(epdata)) stop("epdata must be a data.frame",call.=FALSE)
  if (!all(names(epdata)[1:7]==c("id","Tstart","Tstop","status","agestart",
    "agestop","birth")))
        stop ("epdata must have: Tstart,Tstop,status,agestart,agestop,birth",call.=FALSE) 
  if (!missing(weight)) {if (!is.numeric(weight)) stop ("weight must be numeric",call.=FALSE)}
  if (!missing(formula) & formula!=~1) {
     covar <- model.frame(formula=formula, data=epdata)     
     factors<-names(covar)
     nlev<-0
     for (k in 1:length(factors)) {
        lev.labels<-
        nlev<-nlev+nlevels(covar[,k])
        }
  }
  if (!missing(win)) {if (any(is.na(win))) stop("NA win for one or more obs.",call.=FALSE)}
  if (method!="cpr" & method!="car" & method!="cap")
      stop("Method must be 'cpr','car', or 'cap'",call.=FALSE)
  if (agelimits[1]>=agelimits[2] | agelimits[1]<0 )
      stop("Unconsistent agelimits for one or more obs",call.=FALSE)
  
  # weight
  if (!missing(weight)) epdata$weight<-weight
  if (missing(weight)) epdata$weight<-rep(1,dim(epdata)[1])

  # Calendar Window of observation______________________________________________
  if (!missing(win)) {
     epdata$Wstart<-pmax(win[,1],epdata$Tstart,na.rm=TRUE)
     epdata$Wstop<-pmin(win[,2],epdata$Tstop, na.rm=TRUE)
     chk<-ifelse(epdata$Wstop<epdata$Wstart,0,1)
  } else {
     epdata$Wstart<-epdata$Tstart
     epdata$Wstop<-epdata$Tstop
     chk<-1
  }
  
  epdata<-subset(epdata,subset & chk==1)
  chk<-NULL
  epdata$Astart<-pmax(epdata$Wstart-epdata$birth)
  epdata$Astop<-pmin(epdata$Wstop-epdata$birth)

  # Check number of events
  if (dim(epdata)[1]==0 | all(epdata$status==0)) stop("No events in the dataframe")
  
  # BASE MODEL (NO COVARIATES)__________________________________________________
  require(mgcv, quietly=TRUE)
  # Aggregate data (covariate is constant 1)
  
  a<-aggreg(epdata,method,rep(1,dim(epdata)[1]))
  #select age range and rows with non-zero exposure time
  a<-subset(a,a$age>=agelimits[1] & a$age<=agelimits[2] & a$expos>0.001)
  age0<-a$age
  event0<-a$event
  expos0<-a$expos
  unsmoo<-cbind(age0,round(a$event/a$expos,4))
  colnames(unsmoo)<-c("age","base")
  # GAM Model
  m0<-gam(round(event0)~offset(log(expos0))+ s(age0),family=poisson)
  baseline0<-as.vector(exp(predict.gam(m0,type="terms",terms="s(age0)")+m0$coeff[1]))

  # Check number of events
  if (sum(event0)<=50)
      warning("The number of events in the window of observation is lower than 50.",
              immediate.=TRUE,call.=FALSE)
  # Tail flattening
  if (any(tails)) baseline0<-tailflat(age0,baseline0,cumsum(event0/sum(event0)),tails)
  tab<-cbind(age0,round(baseline0,7))
  colnames(tab)<-c("age","baseline")


  # MODEL WITH COVARIATES_______________________________________________________
  if  (formula!=~1) {
      # Initialize objects
      catcnt<-1
      anovatest<-rep(NA,length(factors))
      # initialize (empty) arrays for final output
      outKNOT<-matrix(NA,0,2)
      outBETA<-matrix(NA,0,4)
      outEVENT<-matrix(NA,0,4)
      outSE<-matrix(NA,0,4)
      outPVALUE<-matrix(NA,0,4)

      # select k-th covariate
   for (k in 1:length(factors)) {    # k-th covariate      
      # Aggregate data                
      ak<-aggreg(epdata,method,epdata[,factors[k]])
      #select age range and rows with non-zero exposure time
      ak<-subset(ak,ak$age>=agelimits[1] & ak$age<=agelimits[2] &
           ak$expos>0.001)
      agek<-ak$age
      eventk<-ak$event
      exposk<-ak$expos
      evexp<-eventk/exposk              # UNSMOOTHED RATES
      variab<-factor(ak$variab)
      ncat<-nlevels(variab)  # number of levels of covariate
      categ<-levels(variab)  # levels of covariate
      if (k==1) lev.labels<-categ
      if  (k>1) lev.labels<-c(lev.labels,categ)
      
      for (z in 1:ncat)    {    # z-th level
           rrisk<-rep(0,length(age0))
           # defining knots and age intervals (knots are at exact ages)
           prop<-cumsum(evexp[variab==categ[z]]/sum(evexp[variab==categ[z]]))
           knot<-c(0,0)
           knot[1]<- min(agek[variab==categ[z]][which(prop>=0.33)])
           knot[2]<- min(agek[variab==categ[z]][which(prop>=0.667)])
           # mid-points between jumps and distances between mid-points
           midpt<-rep(0,3)
           midpt[1]<- min(agek[variab==categ[z]][which(prop>=0.15)])
           midpt[2]<- min(agek[variab==categ[z]][which(prop>=0.50)])
           midpt[3]<- min(agek[variab==categ[z]][which(prop>=0.85)])
           ddist<-c(0,0)
           ddist[1]<-which(agek[variab==categ[z]]==midpt[2])-
                         which(agek[variab==categ[z]]==midpt[1])
           ddist[2]<-which(agek[variab==categ[z]]==midpt[3])-
                         which(agek[variab==categ[z]]==midpt[2])

           int<-as.factor(ifelse(agek<knot[1],1,
                         ifelse(agek>=knot[1] & agek<=knot[2],2,3)))
           # 'deviation coding' dummies (for each level and age subint.)
           dummy<-matrix(rep(0,3*(ncat-1)*dim(ak)[1]),ncol=3*(ncat-1))
           j<-0
           if (z<ncat){
              for (i in 1:3) {           # i-th age subinterval
                for (l in 1:(ncat-1)) { # l-th category
                   j<-j+1               # j-th dummy
                   dummy[,j]<-ifelse(variab==categ[l] & int==i,1,
                              ifelse(variab==categ[ncat] & int==i,-1,0))}}}
              j<-0
           if (z==ncat) {
              for (i in 1:3) {           # i-th age subinterval
                for (l in 1:(ncat-1)) { # l-th category
                   j<-j+1               # j-th dummy
                   dummy[,j]<-ifelse(variab==categ[l+1] & int==i,1,
                              ifelse(variab==categ[1] & int==i,-1,0))}}}


           # Creates a formula for the model
           fmla <- as.formula(paste("round(eventk)~offset(log(exposk))+ s(agek)+",
                paste("dummy[,",1:(3*(ncat-1)),"]",sep="",collapse="+"),sep=""))
           # Model
           mk<-gam(fmla,family=poisson)
           baselinek<-exp(predict.gam(mk,type="terms",terms="s(agek)")+mk$coeff[1])

           #relative risks [row=categories; col=age interval]
           if (z<ncat) {
                betax<-c(exp(mk$coeff[z+1]),
                         exp(mk$coeff[ncat+z]),
                         exp(mk$coeff[2*ncat+z-1]),NA)
                 se<-c(summary.gam(mk)$se[z+1],
                       summary.gam(mk)$se[ncat+z],
                       summary.gam(mk)$se[2*ncat+z-1],NA)
                 pvalue<-c(summary.gam(mk)$p.pv[z+1],
                       summary.gam(mk)$p.pv[ncat+z],
                       summary.gam(mk)$p.pv[2*ncat+z-1],NA)
           }
           if (z==ncat) {
                betax<-c(exp(mk$coeff[ncat]),
                         exp(mk$coeff[2*ncat-1]),
                         exp(mk$coeff[3*ncat-2]),NA)
                se<-c(summary.gam(mk)$se[ncat],
                       summary.gam(mk)$se[2*ncat-1],
                       summary.gam(mk)$se[3*ncat-2],NA)
                pvalue<-c(summary.gam(mk)$p.pv[ncat],
                       summary.gam(mk)$p.pv[2*ncat-1],
                       summary.gam(mk)$p.pv[3*ncat-2],NA)
           }

           # vectors of baseline without duplications
           baseline<-rep(0,length(age0))
           j<-1
           baseline[j]<-baselinek[1]
           for (i in 1:length(agek)) {
              if (agek[i]!=age0[j]) {
                   j<-j+1
                   baseline[j]<-baselinek[i]
              }
           }
           # SMOOTHING PROCEDURE (from step-curves to continous curves)_______
           # using logistic weights
           w1<-c(rep(0,ddist[1]))
           w2<-c(rep(0,ddist[2]))
           if (sig.eff) bc<-ifelse(pvalue<sig.lev,betax,1)
           if (!sig.eff) bc<-betax
           r<-5/ddist
           cc<-exp(((ddist+1)*r)/2)
           w1<-1/(1+cc[1]*exp(-r[1]*seq(1:ddist[1])))
           w2<-1/(1+cc[2]*exp(-r[2]*seq(1:ddist[2])))
           rrisk[age0<midpt[1]]<-bc[1]
           rrisk[age0>=midpt[1] & age0<midpt[2]]<-bc[1]*(1-w1)+bc[2]*w1
           rrisk[age0==midpt[2]]<-bc[2]
           rrisk[age0>midpt[2] & age0<=midpt[3]]<-
                bc[2]*(1-w2)+bc[3]*w2
           rrisk[age0>midpt[3]]<-bc[3]
           rrisk<-rrisk*baseline
           if (any(tails)) rrisk<-tailflat(age0,rrisk,prop,tails)

           # Preparing output
           tab<-cbind(tab,rrisk) # Age profiles
           outKNOT<-rbind(outKNOT,knot)  # knots
           outEVENT<-rbind(outEVENT,c(   # events
             sum(eventk[variab==categ[z] & int==1]),
             sum(eventk[variab==categ[z] & int==2]),
             sum(eventk[variab==categ[z] & int==3]),
             sum(eventk[variab==categ[z]])))           
           outBETA<-rbind(outBETA,betax)   # relative risks
           outSE<-rbind(outSE,se)  # Standard errors
           outPVALUE<-rbind(outPVALUE,pvalue)  # P-values
           temp<-rep(0,length(age0))
           for (i in 1:length(age0)){
               if (any(agek==age0[i] & variab==categ[z]))
                     {temp[i]<-evexp[agek==age0[i] & variab==categ[z]]}}
           unsmoo<-cbind(unsmoo,temp)
           catcnt<-catcnt+1 # Counter
      } # close iteration for z-th level
      
      # Total column (Model without age intervals)
      contrasts(variab)<-'contr.sum'
      m1<-gam(round(eventk)~offset(log(exposk))+ s(agek)+variab,family=poisson)
       
      # changes the order of levels
      last<-factor(as.character(variab), levels=c(categ[ncat],categ[1:(ncat-1)]))
      contrasts(last)<-'contr.sum'
      m2<-gam(round(eventk)~offset(log(exposk))+ s(agek)+last,family=poisson)

      outBETA[(catcnt-ncat):(catcnt-1),4]<-
                c(exp(m1$coeff[2:ncat]),exp(m2$coeff[2]))
      outSE[(catcnt-ncat):(catcnt-1),4]<-
                c(summary.gam(m1)$se[2:ncat],summary.gam(m2)$se[2])
      outPVALUE[(catcnt-ncat):(catcnt-1),4]<-
                c(summary.gam(m1)$p.pv[2:ncat],summary.gam(m2)$p.pv[2])

      # ANOVA test (comparison between base model and model with covariate__
      m0k<-gam(round(eventk)~offset(log(exposk))+ s(agek),family=poisson)
      test<-anova.gam(m0k,m1,test="Chisq")
      anovatest[k]<-round(test$P[2],3)  # ANOVA test pvalue

   }  # close iteration for k-th covariate
   # FINAL OUTPUT
   rownames(outBETA)<-lev.labels
   rownames(outEVENT)<-lev.labels
   rownames(outSE)<-lev.labels
   rownames(outPVALUE)<-lev.labels
   rownames(outKNOT)<-lev.labels
   colnames(tab)<-c("age","baseline",lev.labels)
   colnames(unsmoo)<-c("age","baseline",lev.labels)
   colnames(outKNOT)<-c("knot1","knot2")
   colnames(outBETA)<-c("int1","int2","int3","tot")
   colnames(outEVENT)<-c("int1","int2","int3","tot")
   colnames(outSE)<-c("int1","int2","int3","tot")
   colnames(outPVALUE)<-c("int1","int2","int3","tot")
   outlist<-list(tab,unsmoo,outKNOT, round(outEVENT),round(outBETA,4),
             round(outSE,4),round(outPVALUE,3),factors,anovatest,tr.name,method,
             agelimits,tails)
   names(outlist)<-c("profiles","unsmoothed","knot","event","rrisk","se",
         "pvalue","factors","ANOVAtest","tr.name","method","agelimits","tails")
   }  # close non-missing formula
   else
   {
      outlist<-list(tab,unsmoo,tr.name,method,agelimits,tails)
      names(outlist)<-c("profiles","unsmoothed","tr.name","method",
                        "agelimits","tails")
   }
   # WRITE TEXT FILE
   if (outfile) {
      sink(paste(tr.name,".txt",sep=""))
      print(outlist)
      sink()
   }
   if (any(outEVENT<5)) 
        warning("Few events (<5) for one or more levels of specified factors. Check the output carefully"
            ,call.=FALSE)   
   class(outlist)<-"ap"
   return(outlist)
}




#_______________________________________________________________________________
# epdata   creates episode data structure from dates
#_______________________________________________________________________________

epdata<-function(start, event, lcensor, rcensor,subset=TRUE, birth,id, addvar)
{
  #  CHECK arguments
  if (missing(start))    stop("Must have a start argument.",call.=FALSE)
  if (missing(event))    stop("Must have a event argument.",call.=FALSE)
  if (!is.logical(subset)) stop("subset is not logical.",call.=FALSE)
  if (missing(id))       id<- seq(1,length(start))
  if (missing(lcensor))  lcensor<- rep(NA,length(start))
  if (missing(rcensor))  rcensor<- rep(NA,length(start))
  if (length(start)!=length(event) | length(start)!=length(id))
        stop ("Arguments must have same length.",call.=FALSE)
  if (!missing(addvar)){
    if (dim(as.matrix(addvar))[1]!=length(start))
        stop ("Arguments must have same length.",call.=FALSE)
     addvar<-subset(addvar,subset)
  }
  id<-id[subset]
  if (any(is.na(id))) stop("NA id for one or more obs.",call.=FALSE)
  start<-start[subset]
  if (!is.numeric(start)) stop ("Start variable is not numeric.",call.=FALSE)
  event<-event[subset]
  if (!is.numeric(event)) stop ("Event variable is not numeric.",call.=FALSE)
  birth<-birth[subset]
  if (any(is.na(birth))) stop("NA date of birth for one or more obs.",call.=FALSE)
  lcensor<-lcensor[subset]
  rcensor<-rcensor[subset]
  if (any(is.na(start) & is.na(lcensor)))
         stop("Both start and lcensor are NA for some obs.",call.=FALSE)
  if (any(is.na(event) & is.na(rcensor)))
         stop("Both event and rcensor are NA for some obs.",call.=FALSE)

  # episode (start and stop)
  Tstart<-pmax(lcensor,start,na.rm=TRUE)
  Tstop<-pmin(rcensor,event,na.rm=TRUE)
  # STATUS variable 0:right cens; 1:event; 2:left cens; 3:interval cens.
  status<- ifelse(is.na(event),
            ifelse(is.na(lcensor),0,3),
            ifelse(is.na(lcensor),1,2))
  # Ages
  agestart<-Tstart-birth    # age at the beginning of the episode
  agestop<-Tstop-birth      # age at the end of the episode

  # Final dataframe
  e<- cbind(id,Tstart,Tstop,status,agestart,agestop,birth)
  if (!missing(addvar)) e<-cbind(e,addvar)
  return (e)
}


#_______________________________________________________________________________
# splitter   creates time varying variables through episode splitting
#_______________________________________________________________________________


splitter<-function(epdata,split,tvar.name="Tvar",tvar.lev)
{
  # checks arguments____________________________________________________________
  if (missing(epdata)) stop("Must have a epdata argument.",call.=FALSE)
  if (!is.data.frame(epdata)) stop("epdata must be a data.frame.",call.=FALSE)
if (!all(names(epdata)[1:7]==c("id","Tstart","Tstop","status","agestart",
    "agestop","birth")))
        stop ("epdata must have: Tstart,Tstop,status,agestart,agestop,birth",call.=FALSE)
  if (missing(split))  stop("Must have a split argument.",call.=FALSE)
  N<-ifelse(is.vector(split),1,dim(split)[2])
  len<-ifelse(is.vector(split),length(split),dim(split)[1])
  if (dim(epdata)[1]!=len) stop ("epdata and split have different length",call.=FALSE)

  if (missing(tvar.lev)) tvar.lev<-seq(1,N+1)
  if (length(tvar.lev)!=N+1) stop("Unconsistent tvar.lev",call.=FALSE)
  if (N>1) {
    for (k in 1:(N-1)) {
       if (any(!is.na(split[,k]) & !is.na(split[,k+1]) & split[,k]>split[,k+1]))
            stop("Multiple split dates must be sequentials",call.=FALSE)
       if (any(is.na(split[,k]) & !is.na(split[,k+1])))
            stop("In multiple split dates a split date cannot be followed by non-NA date"
              ,call.=FALSE)
    }
  }
  # Prepare wide format (subepisodes on the same record)________________________
  if (N==1)  epdata[,"tvar"]<-split
  if (N>1) {for (k in 1:(N)) {epdata[,paste("tvar",k,sep="")]<-split[,k]}}
  co<-dim(epdata)[2]
  epdata$Ts1<-epdata$Tstart # Tstart of 1st subepisode
  for (k in 1:N) {
      if (N>1) splitk<-split[,k]
      if (N==1) splitk<-split

      # Tstop of k-th subepisode
      epdata[,paste("Tf",k,sep="")]<-
          ifelse(!is.na(splitk) & splitk<epdata$Tstop & splitk>epdata[,paste("Ts",k,sep="")],
              splitk,ifelse(is.na(splitk) | splitk<epdata[,paste("Ts",k,sep="")],
              epdata[,paste("Ts",k,sep="")],epdata$Tstop ))
      # status of k-th subepisode
      epdata[,paste("st",k,sep="")]<-ifelse(!is.na(splitk) & splitk<epdata$Tstop
             & splitk>epdata[,paste("Ts",k,sep="")],
             ifelse(epdata$status<=1,0,3),epdata$status)
      # Tstart of (k+1)-th subepisode
      epdata[,paste("Ts",k+1,sep="")]<-epdata[,paste("Tf",k,sep="")]
  }
  epdata[,paste("Tf",N+1,sep="")]<-epdata$Tstop  # Tstop of (k+1)-th subepisode
  epdata[,paste("st",N+1,sep="")]<-epdata$status # Status of (k+1)-th subepisode

  # create long format: one record for each subepisode__________________________
  elong<-reshape(epdata, direction="long", varying=names(epdata)[(co+1):(co+3+3*N)],
             sep="",idvar="idt")
  elong<-elong[order(elong$id),] # sort cases by id
  elong<-subset(elong, elong$Ts<elong$Tf) # remove duplicated cases

  # time varying variable_______________________________________________________
  elong[,tvar.name]<-tvar.lev[1]
  if (N==1) {
     elong[,tvar.name]<-ifelse(!is.na(elong$tvar) &
                         elong$Ts>=elong$tvar,tvar.lev[2],elong[,tvar.name])
     elong$tvar<-NULL
  }
  if (N>1){
      for (k in 1:N){
          elong[,tvar.name]<-ifelse(!is.na(elong[,paste("tvar",k,sep="")]) &
                elong$Ts>=elong[,paste("tvar",k,sep="")],tvar.lev[k+1],elong[,tvar.name])
      elong[,paste("tvar",k,sep="")]<-NULL
      }
  }
  # arranging final dataframe___________________________________________________
  elong$time<-NULL
  elong$idt<-NULL
  row.names(elong)<-NULL
  elong$Tstart<-elong$Ts
  elong$Tstop<-elong$Tf
  elong$status<-elong$st
  elong$Ts<-NULL
  elong$Tf<-NULL
  elong$st<-NULL
  elong$agestart<-elong$Tstart-elong$birth
  elong$agestop<-elong$Tstop-elong$birth
  elong[,tvar.name]<-factor(elong[,tvar.name])
  return(elong)
}


#_______________________________________________________________________________
# plotap     plots ageprofile
#_______________________________________________________________________________

plotap<-function(x,lev.labels, baseline=TRUE,unsmoothed=FALSE,
           xlim, ylim, title)
{
  # Check arguments ______________________________________________________________
  if (missing(x)) stop ("Must have an ap argument",call.=FALSE)
  if (class(x)!="ap")
        stop("Must have an ap argument",call.=FALSE)
  if (missing(title)) title<-x$tr.name
  age<-x$profiles[,1]
  base<-x$profiles[,2]
   if (!missing(lev.labels)) {
     levsp<-strsplit(lev.labels,split="*",fixed=TRUE)
     rr<-data.frame(matrix(1,length(age),length(lev.labels)))
     urr<-data.frame(matrix(1,length(age),length(lev.labels)))
     for (k in 1:length(lev.labels)) {
        for (j in 1: length(levsp[[k]])) {
           if (any(levsp[[k]][[j]]==colnames(x$profiles))) {
               rr[,k]<-rr[,k]*x$profiles[,levsp[[k]][[j]]]/base
               urr[,k]<-urr[,k]*x$unsmoothed[,levsp[[k]][[j]]]/base
               }
           else {
               stop(paste(levsp[[k]][[j]],"not found in the dataframe."),call.=FALSE)
           }
        }
     }
  }
  if (missing(ylim)) {ylim<-max(base)*2}
  if (missing(xlim)) {xlim<-c(min(age),max(age))}
  if (x$method=="cpr") ylabel<-"Cohort-period rates"
  if (x$method=="car") ylabel<-"Cohort-age rates"
  if (x$method=="cap") ylabel<-"Cohort-age probs"
  # Open window and plot baseline_______________________________________________
  dev.new()
  plot(c(xlim[1],xlim[2]),c(0,ylim),
          xlab="Age",ylab=ylabel,type="n")
  title(main = title, font.main=3)
  if (baseline) lines(age+0.5,base,lty=1,lwd=2)

  # plot age profiles for each factor___________________________________________
  if (!missing(lev.labels)){
      names(rr)<-lev.labels
      for (k in 1:length(lev.labels)) {
        lines(age+0.5,rr[,k]*base,col=k,lty=k+1)
        if (unsmoothed) {points(age+0.5,urr[,k]*base, col=k)}
      }
      # Plot legend
      posx<-0.75* xlim[2]+0.25*xlim[1]
      ce<-2-((30+max(nchar(lev.labels)))/37)
      legend (posx,ylim,lev.labels,lty=2:(length(lev.labels)+1),col=1:length(lev.labels),
       cex=ce,y.intersp=0.8)
  }
  if (missing(lev.labels) & unsmoothed) points(age+0.5,x$unsmoothed[,2])
}


#_______________________________________________________________
# tabx     Print univariate and bivariate frequency distribution
#_______________________________________________________________     
tabx<-function(x,y,prow=FALSE,pcol=FALSE, chisq=FALSE)
{ 
 # uni-variate table
 if (missing(y)) 
  {
    n<-table(x,useNA="ifany")
    N<-sum(n)
    perc<-round(n/N,2)          
        
    if (any(is.na(x)))   # with NA values
    { 
      valid<-c(round(table(x)/sum(table(x)),2),"")
      Tot<-c(sum(n),1,1)
      tabf<-rbind(cbind(n,perc,valid),Tot)
     } else             # no NA values
    { 
      Tot<-c(sum(n),sum(perc))
      tabf<-rbind(cbind(n,perc),Tot)
     }  
    print("Univariate frequency table",quote=FALSE)
    print(tabf,quote=FALSE)
  }
  # bivariate table
  if (!missing(y)) {
    tabb<-table(x,y,useNA="ifany")
    N<-sum(tabb) 
    # 1. Counts  
    Tot<-margin.table(tabb,1)
    tot.c<-margin.table(tabb,2)
    tabf<-rbind(cbind(tabb,Tot),c(tot.c,N))  
    rownames(tabf)[dim(tabb)[1]+1]<-"Tot" 
    print("Bivariate frequency table",quote=FALSE)
    print("Counts",quote=FALSE)  
    print(tabf,quote=FALSE)
    if (chisq) print(summary(tabb))
    if (prow) {                   # row percentages      
      n<-margin.table(tabb,1)
      tab1<-prop.table(tabb,1)
      Tot<-margin.table(tab1,1)
      tot.c<-round(margin.table(tabb,2)/N,2)
      tabf<-rbind(cbind(tab1,Tot,n),c(tot.c,sum(tot.c),N))  
      rownames(tabf)[dim(tabb)[1]+1]<-"Tot"
      print("Row percentages",quote=FALSE)  
      print(round(tabf,2),quote=FALSE)
    } 
    if (pcol) {                   # column percentages
      n<-margin.table(tabb,2)
      tab1<-prop.table(tabb,2)
      Tot<-margin.table(tab1,2)
      tot.r<-margin.table(tabb,1)/N
      tabf<-cbind(rbind(tab1,Tot,n),c(tot.r,sum(tot.r),N))  
      colnames(tabf)[dim(tabb)[2]+1]<-"Tot"
      print("Column percentages",quote=FALSE)
      print(round(tabf,2),quote=FALSE) 
    }      
  }

}

#_______________________________________________________________
# tabm     Print multiple regression and logit model estimations
#_______________________________________________________________                 
tabm<-function(mod,pvalues=TRUE,digits=3)
{
  Beta<-round(coef(mod),digits)
  sign1<-summary(mod)$coef[,4]
  sign2<-as.character(cut(sign1,
    breaks=c(0,0.01,0.05,0.1,1),
    labels=c("***","**","*","")))    
  if (attributes(mod)$class[1]=="glm"){
      model<-"Logistic regression."
      fit<-paste("AIC:",round(AIC(mod),2))
      expBeta<-as.character(round(exp(Beta),digits))
      if (names(coef(mod))[1]=="(Intercept)") expBeta[1]<-"--" 
        if (pvalues){
          tabb<-cbind(Beta,expBeta,round(sign1,digits),sign2) 
          colnames(tabb)<-c("Beta","exp(Beta)","sign.","")
        } else {
          tabb<-cbind(Beta,expBeta,sign2)
          colnames(tabb)<-c("Beta","exp(Beta)","sign.")
        }
  }
  if (attributes(mod)$class[1]=="lm"){
      model<-"Multiple regression."
      fit<-paste("AIC:",round(AIC(mod),2),"   Adjusted R square:",
            round(summary(mod)$adj.r.squared,2))
      if (pvalues){  
          tabb<-cbind(Beta,round(sign1,digits),sign2) 
          colnames(tabb)<-c("Beta","sign.","") 
      } else {             
          tabb<-cbind(Beta,sign2)
          colnames(tabb)<-c("Beta","sign.")
      }
  }
  print (model,quote=FALSE)
  print(tabb,quote=FALSE)
  print("Sign.:'***' <0.01; '**' <0.05; '*' <0.1",quote=FALSE)
  print(fit,quote=FALSE)
}


                 
#_______________________________________________________________________________
# mkdate: creates dates in continuous years or cmc
#_______________________________________________________________________________

# Dates of event must be (year: XXXX format; month XX format);
# option: cmc=FALSE years continuous time - cmc=TRUE date in Century month code)

mkdate<-function(year,month,cmc=FALSE)
{
  if (!cmc)   mkdate<-year-1900+round((month-0.5)/12,2)
  if (cmc)    mkdate<-(year-1900)*12+month
  return(mkdate)
}

                
#_______________________________________________________________________________
# listvar: # Show variables in a dataframe and the associate number of columns
#_______________________________________________________________________________

listvar=function(df)
{      
    temp1=names(df)
    temp2=seq(1,length(temp1))
    temp3=data.frame(temp1,temp2)
    names(temp3)=c("VAR","COL")
    return(temp3)
    rm(temp1,temp2,temp3)
}



  