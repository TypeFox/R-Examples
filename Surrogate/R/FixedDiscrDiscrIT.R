
FixedDiscrDiscrIT<-function(Dataset,Surr,True,Treat,Trial.ID,Weighted=TRUE,Setting=c("binord")){
  
#	if (!requireNamespace("OrdinalLogisticBiplot", quietly = TRUE)|!requireNamespace("logistf", quietly = TRUE)) {
#    stop("OrdinalLogisticBiplot and logistf packages needed for this function to work. Please install them.",
#      call. = FALSE)  }

#	if (!requireNamespace("rms", quietly = TRUE)|!requireNamespace("MASS", quietly = TRUE)) {
#	  stop("rms and MASS packages needed for this function to work. Please install them.\r\n",
#	       call. = FALSE)  }

#library(OrdinalLogisticBiplot)
#library(logistf)

#library(rms)
#library(MASS)

# removing missing data
Total.length<-length(Dataset[,1])  
Dataset<-na.omit(Dataset)
Removed.obs<-Total.length-length(Dataset[,1])
if (Removed.obs>1) {cat("Problem: parewise deletion due to missing data in", Removed.obs,"rows")}


  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat<- Dataset[,paste(substitute(Treat))]
  Trial.ID <- Dataset[,paste(substitute(Trial.ID))]

if (Setting==c("binord")|Setting==c("ordord")) {
while (any(True<=0)) {True<-True+1}
}
if (Setting==c("ordbin")|Setting==c("ordord")) {
while (any(Surr<=0)) {Surr<-Surr+1}
}
  g2<-NA;l.ncpara<-NA;u.ncpara<-NA;Obs.per.trial<-Trial<-NULL;r2hmax<-NA;l.r2h<-NA;u.rhmax<-NA;expg2<-0;error1<-error2<-NA
  


trialnumber<-length(unique(Trial.ID))
R2h.included<-rep("Y",trialnumber)
R2ht.included<-rep("Y",trialnumber)

  ######################################################
  ###### R2h.max Individual level surrogacy   ##########
  ######################################################
  
  
  if (Setting==c("binord")){
    for (i in 1:trialnumber){
      Obs.per.trial[i]<-length(Treat[Trial.ID==i])
      Trial[i]<-i
      # R2H CODE lrm within loop
      error1[i]<-0
      error2[i]<-0
      tryCatch({mod1<-rms::lrm(True[Trial.ID==i]~Treat[Trial.ID==i])},error=function(e){error1[i]<<-1},warning=function(w){})
      tryCatch({mod2<-rms::lrm(True[Trial.ID==i]~Treat[Trial.ID==i]+Surr[Trial.ID==i])},error=function(e){error2[i]<<-1},warning=function(w){})
      if (error1[i]!=1&error2[i]!=1){
        g2[i]<-2*(as.numeric(logLik(mod2))-as.numeric(logLik(mod1)))
        expg2<-expg2+exp(-g2[i]/Obs.per.trial[i])
        
        #### R2H CODE lrm confidence intervals
        if (is.nan(g2[i])) {g2[i]<-NA; {cat("Problem R2h: model failure for one of the models of trial", i, "/N.\r\n")}
		R2h.included[i]<-"N"} else {
          if (abs(g2[i])<1e-04) {g2[i]<-0}
          l.ncpara[i]<-qchisq(0.025,1,g2[i])
          if (pchisq(g2[i],1,0)<=0.95) {l.ncpara[i]<-0}
          u.ncpara[i]<-qchisq(0.975,1,g2[i])
        }}else{
          if (error1[i]!=1){
            # if there are no surrogate outcomes then mod2 will fail to run in this case g2 is set to 0
            g2[i]<-0
            expg2<-expg2+exp(-g2[i]/Obs.per.trial[i])
            
            #### R2H CODE lrm confidence intervals
            if (is.nan(g2[i])) {g2[i]<-NA; {cat("Problem R2h: model failure for one of the models of trial", i, "/N.\r\n")}
			R2h.included[i]<-"N"} else {
              # sometimes difference is negative but so close to zero as to make no odds to allow these trials to contibute to CIs recoded to be 0
              if (abs(g2[i])<1e-04) {g2[i]<-0}
              l.ncpara[i]<-qchisq(0.025,1,g2[i])
              if (pchisq(g2[i],1,0)<=0.95) {l.ncpara[i]<-0}
              u.ncpara[i]<-qchisq(0.975,1,g2[i])}}else {{cat("Problem R2h: model failure for one of the models of trial", i, "/N.\r\n")}
			R2h.included[i]<-"N"}
        }}
    modRmax<-rms::lrm(True~1)
  }
  if (Setting==c("ordbin")){
    
    for (i in 1:trialnumber){
      # R2H CODE lrm within loop
	error1[i]<-0
      error2[i]<-0
      Obs.per.trial[i]<-length(Treat[Trial.ID==i])
      Trial[i]<-i
      tryCatch({mod1<-glm(True[Trial.ID==i]~Treat[Trial.ID==i],family=binomial(link=logit))},error=function(e){error1[i]<<-1})
      tryCatch({mod2<-glm(True[Trial.ID==i]~Treat[Trial.ID==i]+Surr[Trial.ID==i],family=binomial(link=logit))},error=function(e){error2[i]<<-1})
      if (error1[i]!=1&error2[i]!=1){
        g2[i]<-2*(as.numeric(logLik(mod2))-as.numeric(logLik(mod1)))
        expg2<-expg2+exp(-g2[i]/Obs.per.trial[i])
        
        # R2H CODE lrm confidence intervals
        if (is.nan(g2[i])) {g2[i]<-NA; {cat("Problem R2h: model failure for one of the models of trial", i, "/N.\r\n")}
		R2h.included[i]<-"N"} else {
          l.ncpara[i]<-qchisq(0.025,1,g2[i])
          if (pchisq(g2[i],1,0)<=0.95) {l.ncpara[i]<-0}
          u.ncpara[i]<-qchisq(0.975,1,g2[i])
        }}else {{cat("Problem R2h: model failure for one of the models of trial", i, "/N.\r\n")}
		R2h.included[i]<-"N"}
    }
    modRmax<-glm(True~1,family=binomial(link=logit))
  }
  if (Setting==c("ordord")){
    
    for (i in 1:trialnumber){
      # R2H CODE lrm within loop
      error1[i]<-0
      error2[i]<-0
      Obs.per.trial[i]<-length(Treat[Trial.ID==i])
      Trial[i]<-i
      tryCatch({mod1<-rms::lrm(True[Trial.ID==i]~Treat[Trial.ID==i])},error=function(e){error1[i]<<-1},warning=function(w){})
      tryCatch({mod2<-rms::lrm(True[Trial.ID==i]~Treat[Trial.ID==i]+as.numeric(Surr[Trial.ID==i]))},error=function(e){error2[i]<<-1},warning=function(w){})
      if (error1[i]!=1&error2[i]!=1){
        g2[i]<-2*(as.numeric(logLik(mod2))-as.numeric(logLik(mod1)))
        expg2<-expg2+exp(-g2[i]/Obs.per.trial[i])
        
        #### R2H CODE lrm confidence intervals
        if (is.nan(g2[i])) {g2[i]<-NA; {cat("Problem R2h: model failure for one of the models of trial", i, "/N. \r\n ")}
		R2h.included[i]<-"N"} else {
          if (abs(g2[i])<1e-04) {g2[i]<-0}
          l.ncpara[i]<-qchisq(0.025,1,g2[i])
          if (pchisq(g2[i],1,0)<=0.95) {l.ncpara[i]<-0}
          u.ncpara[i]<-qchisq(0.975,1,g2[i])
        }}else{{cat("Problem R2h: model failure for one of the models of trial", i, "/N. \r\n ")}
		R2h.included[i]<-"N"}}
    
    modRmax<-rms::lrm(True~1)
    
  }
  
  ## R2Hmax CODE lrm out of loop
  
  lrf<-1-(expg2/length(g2[!is.na(g2)]))
  
  R2h.max<-lrf/(1-exp(2*as.numeric(logLik(modRmax))/sum(Obs.per.trial[R2h.included=="Y"])))
  
  ##R2hmax confidence intervals
  
  paral<-exp(-l.ncpara[!is.na(l.ncpara)]/Obs.per.trial[!is.na(l.ncpara)])
  parau<-exp(-u.ncpara[!is.na(u.ncpara)]/Obs.per.trial[!is.na(u.ncpara)])
  
  lower<-1-(sum(paral)/length(paral))
  upper<-1-(sum(parau)/length(parau))
  
  R2h.max.lb<-lower/(1-exp(2*as.numeric(logLik(modRmax))/sum(Obs.per.trial[R2h.included=="Y"])))
  R2h.max.ub<-min(1, upper/(1-exp(2*as.numeric(logLik(modRmax))/sum(Obs.per.trial[R2h.included=="Y"]))))
  
  
  R2h<-data.frame(cbind(R2h.max, R2h.max.lb, R2h.max.ub)) 
  
  
  
  ############################################################
  ############# R2ht trial level surrogacy ###################
  ############################################################
  
  # record and deal with issues of Separation
level<-unique(Treat)
  Intercept.S<-Treatment.T<-Intercept.T<-Treatment.S<-Separation.S<-Separation.T<-NA
  f.True<-as.numeric(True)
  
  if (Setting==c("binord")){
    for (i in 1:trialnumber){
      FIRTH<-FALSE
      for (S in 0:1){
        for (t in 1:2){
          if (length(Trial.ID[Trial.ID==i&Surr==S&Treat==level[t]])==0){
            # Separation exists firth adjustment required
            FIRTH<-TRUE
          }}}
      
      if (FIRTH==TRUE){
        Separation.S[i]<-"YES"
        modr2ht.firth1<-logistf::logistf(Surr[Trial.ID==i]~Treat[Trial.ID==i],family=binomial(link=logit),pl=FALSE,firth=FIRTH)
        coef.firth<-coefficients(modr2ht.firth1)
        Intercept.S[i]<-coef.firth[1]
        Treatment.S[i]<-coef.firth[2]}else{
          Separation.S[i]<-"NO"
          tryCatch({modr2ht<-glm(Surr[Trial.ID==i]~Treat[Trial.ID==i],family=binomial(link=logit))},error=function(e){})
          coef<-coefficients(modr2ht)
          Intercept.S[i]<-coef[1]
          Treatment.S[i]<-coef[2]}
    }}
  
  
  if (Setting==c("ordbin")){
    
    for (i in 1:trialnumber){
      FIRTH<-FALSE
      for (S in 0:1){
        for (t in 1:2){
          if (length(Trial.ID[Trial.ID==i&True==S&Treat==level[t]])==0){
            # Separation exists firth adjustment required
            FIRTH<-TRUE
          }}}
      if (FIRTH==TRUE){
        Separation.T[i]<-"YES"
        modr2ht.firth1<-logistf::logistf(True[Trial.ID==i]~Treat[Trial.ID==i],family=binomial(link=logit),pl=FALSE,firth=FIRTH)
        coef.firth<-coefficients(modr2ht.firth1)
        Treatment.T[i]<-coef.firth[2]}else{
          Separation.T[i]<-"NO"
          tryCatch({modr2ht.firth1<-glm(True[Trial.ID==i]~Treat[Trial.ID==i],family=binomial(link=logit))},error=function(e){})
          coef.firth<-coefficients(modr2ht.firth1)
          Treatment.T[i]<-coef.firth[2]}}
  }
  
  if (Setting==c("binord")|Setting==c("ordord")){
    Level.True<-unique(True)
    for (i in 1:trialnumber){
      #### TRUE
      modht1.f<-NULL
      ordFIRTH<-FALSE
      
      tryCatch({modht1.f<-polr(True[Trial.ID==i]~Treat[Trial.ID==i],method=("logistic"))},error=function(e){},warning=function(w){})
      if (length(coefficients(modht1.f))==0|is.null(modht1.f)) {Treatment.T[i]<-NA; ordFIRTH<-TRUE}else{
        Treatment.T[i]<-coefficients(modht1.f)}
      
      modr2ht.firth2<-NULL
      
      # recording Separation for ordinal variable
      for (T in 1:length(Level.True)){
        for (t in 1:2){
          # Only one outcome for one treatment in a trial-complete Separation
          if (length(Trial.ID[Trial.ID==i&True==Level.True[T]&Treat==level[t]])!=0&length(Trial.ID[Trial.ID==i&True!=Level.True[T]&Treat==level[t]])==0){
            ordFIRTH<-TRUE
          }}}
      
      #If there are no patients in one treatment group for a trial the code will not work
      if (length(True[Treat == level[1] & Trial.ID == i])==0|length(True[Treat == level[2] & Trial.ID == i])==0){
        cat("Problem: R2ht, only one treatment group in trial",i,"/N, trial removed \r\n")
	R2ht.included[i]<-"N"
        Separation.T[i]<-"NA"
      }else{
        #all the results of one are higher than the other treatment
        if (max(as.numeric(True[Treat==level[1] & Trial.ID==i]))<=min(as.numeric(True[Treat==level[2] & Trial.ID==i]))|max(as.numeric(True[Treat==level[2]& Trial.ID==i]))<=min(as.numeric(True[Treat==level[1]& Trial.ID==i]))){
          ordFIRTH<-TRUE}
        #only run firth if Separation occurs
        if (ordFIRTH==TRUE) {
          Separation.T[i]<-"YES"
          #Pordlogist cannot handle zero categories between treatments that can happen with Separation
          #so recode
          if (max(as.numeric(True[Treat==level[1]&Trial.ID==i]))<min(as.numeric(True[Treat==level[2]&Trial.ID==i]))){
            del<-abs(min(as.numeric(True[Treat==level[2]&Trial.ID==i]))-max(as.numeric(True[Treat==level[1]&Trial.ID==i])))-1
            d<-as.numeric(f.True)
            d[Treat==level[2]&Trial.ID==i]<-as.numeric(d[Treat==level[2]&Trial.ID==i])-del
            f.True[Trial.ID==i]<-as.factor(d[Trial.ID==i])
          }
          if (max(as.numeric(True[Treat==level[2]&Trial.ID==i]))<min(as.numeric(True[Treat==level[1]&Trial.ID==i]))){
            del<-abs(min(as.numeric(True[Treat==level[1] & Trial.ID==i]))-max(as.numeric(True[Treat==level[2] & Trial.ID==i])))-1
            d<-as.numeric(f.True)
            d[Treat==level[1] & Trial.ID==i]<-as.numeric(d[Treat==level[1] & Trial.ID==i])-del
            f.True[Trial.ID==i]<-as.factor(d[Trial.ID==i])
          }
          
          
          #if the min value of True in a trial is >2 pord cannot handle this and we recode   
          if (min(as.numeric(True[Trial.ID==i]))>2) {
            f.True[Trial.ID==i]<-as.factor(as.numeric(f.True[Trial.ID==i])-(as.numeric(min(f.True[Trial.ID==i]))-1))
          }
          for (u in 1:(round(as.numeric(length(unique(True))/2,0)))){
            for (dd in 2:length(unique(True))){
              if (length(f.True[Trial.ID==i&True==dd])==0&length(f.True[Trial.ID==i&f.True<dd])!=0&length(f.True[Trial.ID==i&f.True>dd])!=0){
                g<-as.numeric(f.True)
                g[Trial.ID==i&f.True>dd]<-as.numeric(f.True[Trial.ID==i&f.True>dd])-1
                f.True[Trial.ID==i]<-as.factor(g[Trial.ID==i])
              }}
          }
          modr2ht.firth2<-NULL
          tryCatch({modr2ht.firth2<-OrdinalLogisticBiplot::pordlogist(as.numeric(f.True[Trial.ID==i]),Treat[Trial.ID==i],show = FALSE)},error=function(e){})
          if (is.null(modr2ht.firth2)){Treatment.T[i]<-NA 
{cat("Problem R2ht: Firth model failure for one of the models of trial", i, "/N, trial removed from analysis.\r\n")}
	R2ht.included[i]<-"N"}else{Treatment.T[i]<-coefficients(modr2ht.firth2)
}}else{Separation.T[i]<-"NO"}}
    }}


if (Setting==c("ordbin")|Setting==c("ordord")){
Level.Surr<-unique(Surr)
  # When Firth is applied the proportional odds model produces k-1 surrogate intercept variables (k=levels of ordinal surrogate) 
  # inisialising intercept variables that will be averaged for use at stage two of R2ht
  mui.f.o <- array(NA, c(length(unique(Surr))-1,trialnumber),dimnames = list(C = paste0("mui.f.",1:(length(unique(Surr))-1))))
  f.Surr<-Surr
  for (i in 1:trialnumber){
    ordFIRTH<-FALSE
    modht1.f<-NULL
    tryCatch({modht1.f<-polr(as.factor(Surr[Trial.ID==i])~Treat[Trial.ID==i],method=("logistic"))},error=function(e){})
    if (length(coefficients(modht1.f))==0|is.null(modht1.f)) {
      Treatment.S[i]<- mui.f.o[1:(length(unique(Surr))-1),i]<-NA
      ordFIRTH<-TRUE
    }else{
      Treatment.S[i]<-coefficients(modht1.f)
      mui.f.o[1:(length(unique(Surr))-1),i]<-modht1.f$zeta[1:(length(unique(Surr))-1)]
    }
    
    modr2ht.firth2<-NULL
    
    # recording Separation for ordinal variable
    for (T in 1:length(Level.Surr)){
      for (t in 1:2){
        if (length(Trial.ID[Trial.ID==i&Surr==Level.Surr[T]&Treat==level[t]])!=0&length(Trial.ID[Trial.ID==i&Surr!=Level.Surr[T]&Treat==level[t]])==0){
          Treatment.S[i]<-mui.f.o[1:(length(unique(Surr))-1),i]<-NA
          ordFIRTH<-TRUE
        }}}
    if (length(Surr[Treat==level[1]&Trial.ID==i])==0|length(Surr[Treat==level[2]&Trial.ID==i])==0){
      Treatment.S[i]<-mui.f.o[1:(length(unique(Surr))-1),i]<-NA
      Separation.S[i]<-"NA"
	R2ht.included[i]<-"N"
    } else{
      if (max(as.numeric(Surr[Treat==level[1]&Trial.ID==i]))<=min(as.numeric(Surr[Treat==level[2]&Trial.ID==i]))|max(as.numeric(Surr[Treat==level[2]&Trial.ID==i]))<=min(as.numeric(Surr[Treat==level[1]&Trial.ID==i]))){
        Treatment.S[i]<-mui.f.o[1:(length(unique(Surr))-1),i]<-NA
        ordFIRTH<-TRUE
      }
      # Apply Firth for cases of Separation

      if (ordFIRTH==TRUE){
        Separation.S[i]<-"YES"
        #Pordlogist cannot handle zero categories between treatments in a trial so recode
        if (max(as.numeric(Surr[Treat==level[1]&Trial.ID==i]))<min(as.numeric(Surr[Treat==level[2]&Trial.ID==i]))){
          del<-abs(min(as.numeric(Surr[Treat==level[2]&Trial.ID==i]))-max(as.numeric(Surr[Treat==level[1]&Trial.ID==i])))-1
          d<-as.numeric(f.Surr)
          d[Treat==level[2]&Trial.ID==i]<-as.numeric(d[Treat==level[2]&Trial.ID==i])-del
          f.Surr[Trial.ID==i]<-as.factor(d[Trial.ID==i])
        }
        if (max(as.numeric(Surr[Treat==level[2]&Trial.ID==i]))<min(as.numeric(Surr[Treat==level[1]&Trial.ID==i]))){
          del<-abs(min(as.numeric(Surr[Treat==level[1]&Trial.ID==i]))-max(as.numeric(Surr[Treat==level[2]&Trial.ID==i])))-1
          d<-as.numeric(f.Surr)
          d[Treat==level[1]&Trial.ID==i]<-as.numeric(d[Treat==level[1]&Trial.ID==i])-del
          f.Surr[Trial.ID==i]<-as.factor(d[Trial.ID==i])
        }
        
        #if the min value of Surr in a trial is >2 pord cannot handle this and we recode
        
        if (min(as.numeric(Surr[Trial.ID==i]))>2) {
          f.Surr[Trial.ID==i]<-as.factor(as.numeric(f.Surr[Trial.ID==i])-(as.numeric(min(as.numeric(f.Surr[Trial.ID==i])))-1))
        }
        for (u in 1:(round(as.numeric(length(unique(Surr))/2,0)))){
          for (dd in 2:length(unique(Surr))){
            if (length(f.Surr[Trial.ID==i&Surr==dd])==0&length(f.Surr[Trial.ID==i&f.Surr<dd])!=0&length(f.Surr[Trial.ID==i&f.Surr>dd])!=0){
              g<-as.numeric(f.Surr)
              g[Trial.ID==i&f.Surr>dd]<-as.numeric(f.Surr[Trial.ID==i&f.Surr>dd])-1
              f.Surr[Trial.ID==i]<-as.factor(g[Trial.ID==i])
            }}
        }
        
        tryCatch({modr2ht.firth2<-OrdinalLogisticBiplot::pordlogist(as.numeric(f.Surr[Trial.ID==i]),Treat[Trial.ID==i],show = FALSE)},error=function(e){})
        if (is.null(modr2ht.firth2)){Treatment.S[i]<-mui.f.o[1:(length(unique(Surr))-1),i]<-NA
                                     cat("Problem R2ht: Firth model failure for one of the models of trial", i, "/N, trial removed from analysis.\r\n")
        R2ht.included[i]<-"N"}else{Treatment.S[i]<-coefficients(modr2ht.firth2)
              mui.f.o[1:(length(unique(Surr))-1),i]<-modr2ht.firth2$thresholds[1:(length(unique(Surr))-1)]
              
        }}else{Separation.S[i]<-"NO"}}
    }
  # Mean surrogate intercept value calculated for each trial 
  Intercept.S<-apply(mui.f.o,2,mean,na.rm=TRUE)
}

######################
#### R2ht Stage 2 ####
######################


Trial.Spec.Results <- data.frame(Trial, Obs.per.trial,R2h.included,R2ht.included,Separation.S,Separation.T,Intercept.S,Treatment.S,Treatment.T)
colnames(Trial.Spec.Results) <- c(NULL, "Trial", "Obs.per.trial","R2h.Included","R2ht.Included","Sepatation.S","Separation.T", "Intercept.S", "Treatment.S", "Treatment.T")
rownames(Trial.Spec.Results) <- NULL

lrfht1.f<-NULL
lrfht2.f<-NULL

if (Weighted==TRUE){
  lrfht1.f<-lm(Treatment.T~1,weights=Obs.per.trial)
  lrfht2.f<-lm(Treatment.T~Intercept.S+Treatment.S,weights=Obs.per.trial)
}
if (Weighted==FALSE){
  lrfht1.f<-lm(Treatment.T~1,)
  lrfht2.f<-lm(Treatment.T~Intercept.S+Treatment.S,)
}

remain.trials.f<-min(length(Treatment.T[!is.na(Treatment.T)]),length(Intercept.S[!is.na(Intercept.S)]),length(Treatment.S[!is.na(Treatment.S)]))

# check there are at least 3 trials contributing to analysis
if (remain.trials.f<3|logLik(lrfht2.f)==Inf){r2ht.f<-u.r2ht.f<-l.r2ht.f<-NA
                                             cat("Problem: Not enough trials available for stage two R2ht model \r\n")}else{
                                               
                                               g2ht.f<-2*(as.numeric(logLik(lrfht2.f))-as.numeric(logLik(lrfht1.f)))
                                               
                                               R2ht.f<-1-exp(-g2ht.f/remain.trials.f)
                                               
                                               uhtncp.f<-qchisq(0.975,1,g2ht.f)
                                               lhtncp.f<-qchisq(0.025,1,g2ht.f)
                                               if (pchisq(g2ht.f,1,0)<=0.95) {lhtncp.f<-0}
                                               u.r2ht.f<-min(1,1-exp(-uhtncp.f/remain.trials.f))
                                               l.r2ht.f<-1-exp(-lhtncp.f/remain.trials.f)
                                             }


R2ht<-data.frame(cbind(R2ht.f,l.r2ht.f,u.r2ht.f))

fit <- 
  list(Trial.Spec.Results=Trial.Spec.Results,R2ht=R2ht, R2h=R2h)   

class(fit) <- "FixedDiscrDiscrIT"
fit 

}
