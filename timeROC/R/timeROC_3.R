# {{{ input description :

# T              : vector of observed failure times
# delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
# marker         : vector ofmarker values
# other_markers  : (default is NULL, should be a matrix) other markers that can be associated with the censoring mechanism
# cause          : the value that indicates the main event of interest
# weighting      : (default is "marginal") weighting technique for IPCW : "marginal" for Kaplan-Meier, "cox" for proportional hazards cox model, "aalen" for additive aalen model
# times          : vector of times you want to compute the time dependent AUC.
# ROC            : if TRUE, then save True Positive fraction (Sensitivity) and False Positive fraction (1-Specificity)
#                  for all time in vetor times
#iid             : TRUE or FALSE, indicates if we want to compute the iid representation of the AUC estimator

# }}}
# {{{ output description :

#  1) with competing risks : an object of class "ipcwcompetingrisksROC"
# output of interest #
#TP                  : if ROC=TRUE, the matrix of True Positive fraction (sensitivty),
#                      with  columns correspond to time of vector (ordered) times
#FP_1                : if ROC=TRUE, the matrix of False Positive fraction (1-Specificity), 
#                      where a "control" is defined as an event-free subject at time t,
#                      with  columns correspond to time of vector (ordered) times
#FP_2                : if ROC=TRUE, the matrix of False Positive fraction (1-Specificity), 
#                      where a "control" is defined as any subject that is not a case,
#                      with  columns correspond to time of vector (ordered) times
#AUC_1               : vector of AUC for each times in vector (ordered) times computed with FP_1
#AUC_2               : vector of AUC for each times in vector (ordered) times computed with FP_2
# output for more information #
#times               : the input vector (ordered) times. If there is only one time in the input, 
#                      then 0 is added
#survProb            : the vector of probabilities to be event free (all causes) at each time in 
#                      vector (ordered) times (denominator of False positive fraction) 
#CumulativeIncidence : the vector of the cumulative incidences of the type of event of interest at each time in 
#                      vector (ordered) times
# n                  : number of observations 
# Stats              : matrix with numbers of observed cases and controls (both definitions), and censored before time point,
#                      for each  time in vector (ordered) times
#weights             : an object of class "ipcw" (see "pec" package) with all information about the weights
#computation_time    : the computation time
#iid                 : TRUE or FALSE, indicates if we have computed the iid representation of the AUC estimator
# 2) without competing risks : an object of class "ipcwsurvivalROC"
# output of interest #
#TP                  : if ROC=TRUE, the matrix of True Positive fraction (sensitivty),
#                      with  columns correspond to time of vector (ordered) times
#FP                  : if ROC=TRUE, the matrix of False Positive fraction (1-Specificity), 
#                      where a "control" is defined as an event-free subject at time t,
#                      with  columns correspond to time of vector (ordered) times
#AUC                 : vector of AUC for each times in vector (ordered) times 
#iid                 : TRUE or FALSE, indicates if we have computed the iid representation of the AUC estimator
# output for more information #
#times               : the input vector (ordered) times. If there is only one time in the input, 
#                      then 0 is added
#survProb            : the vector of probabilities to be event free (all causes) at each time in 
#                      vector (ordered) times (denominator of False positive fraction) 
#CumulativeIncidence : the vector of the cumulative incidence of the  event at each time in 
#                      vector (ordered) times
# n                  : number of observations 
# Stats              : matrix with numbers of observed cases and controls, and censored before time point,
#                      for each  time in vector (ordered) times
#weights             : an object of class "ipcw" (see "pec" package) with all information about the weights
#computation_time    : the computation time 

# }}}


timeROC<-function(T,delta,marker,other_markers=NULL,cause,weighting="marginal",times,ROC=TRUE,iid=FALSE){
  # {{{ check some inputs
  if (length(delta)!=length(T) | length(marker)!=length(T) | length(delta)!=length(T)){
    stop("lengths of vector T, delta and marker have to be equal\n") }
  if (missing(times)){
    stop("Choose at least one time for computing the time-dependent AUC\n") } 
  if (!weighting %in% c("marginal","cox","aalen")){
    stop("the weighting argument must be marginal (default), cox or aalen.\n") }  
  if (weighting %in% c("cox","aalen") & !missing(other_markers) & !class(other_markers)=="matrix"){
    stop("argument other_markers must be a matrix\n") }
  if (weighting %in% c("cox","aalen") & !missing(other_markers)){
    if(!nrow(other_markers)==length(marker))  stop("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")
  }
  # }}}
  # {{{ check if there are missing values, and delete rows with missing values
  if (weighting %in% c("cox","aalen") & !missing(other_markers) ){
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker,other_markers)),1,prod))
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
    other_markers<-as.matrix(other_markers[is_not_na,])
  }else{
    is_not_na<-as.logical(apply(!is.na(cbind(T,delta,marker)),1,prod)) 
    T<-T[is_not_na]
    delta<-delta[is_not_na]
    marker<-marker[is_not_na]
  }
  # }}} 
  start_computation_time<-Sys.time()
  # {{{ create some usefull objects
  n<-length(T)
  n_marker<-length(unique(marker))
  n_times<-length(times)
  if (n_times==1){times<-c(0,times)
                  n_times<-2}           # trick to use ipcw.cox() even if there is only one time
  times<-times[order(times)]
  times_names<-paste("t=",times,sep="")
  # }}}
  # {{{ output initialisation
  AUC_1<-rep(NA,n_times)
  AUC_2<-rep(NA,n_times)
  CumInci<-rep(NA,n_times)
  surv<-rep(NA,n_times)
  names(AUC_1)<-times_names
  names(AUC_2)<-times_names
  names(CumInci)<-times_names
  names(surv)<-times_names
  Stats<-matrix(NA,nrow=n_times,ncol=4)
  colnames(Stats)<-c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats)<-times_names
  # }}}
  # {{{  computation of weights (1/2)
  # we need to order to use the pec::ipcw() fonction
  order_T<-order(T)
  T <- T[order_T]
  delta <- delta[order_T]
  marker<- marker[order_T]
  # use ipcw function from pec package
  if(weighting=="marginal"){
    weights <- pec::ipcw(Surv(failure_time,status)~1,data=data.frame(failure_time=T,status=as.numeric(delta!=0)),method="marginal",times=times,subjectTimes=T,subjectTimesLag=1)
  }
  if(weighting=="cox"){
    if (missing(other_markers)){marker_censoring<-marker } 
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="cox",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  if(weighting=="aalen"){
    if (missing(other_markers)){marker_censoring<-marker }
    other_markers<-other_markers[order_T,]
    marker_censoring<-cbind(marker,other_markers)
    colnames(marker_censoring)<-paste("X", 1:ncol(marker_censoring), sep="")
    fmla <- as.formula(paste("Surv(T,status)~", paste(paste("X", 1:ncol(marker_censoring), sep=""), collapse= "+")))
    data_weight<-as.data.frame(cbind(data.frame(T=T,status=as.numeric(delta!=0)),marker_censoring))
    weights <- pec::ipcw(fmla,data=data_weight,method="aalen",times=as.matrix(times),subjectTimes=data_weight[,"T"],subjectTimesLag=1)
  }
  # we order by marker values (in order to compute Se and Sp)
  order_marker<-order(-marker)
  Mat_data<-cbind(T,delta,marker)[order_marker,]
  colnames(Mat_data)<-c("T","delta","marker")
  # Create some weights
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  Weights_cases_all<-Weights_cases_all[order_marker]
  # }}}
  # {{{ Make TP and FP outputs if needed
  if(ROC==TRUE){ 
    FP_1<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    TP<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    FP_2<-matrix(NA,nrow=(n_marker+1),ncol=n_times)
    colnames(FP_1)<-times_names
    colnames(TP)<-times_names
    colnames(FP_2)<-times_names
  } else{FP_1<-NA
         FP_2<-NA
         TP<-NA}
  # }}}
  # {{{ loop on all timepoints t
  for(t in 1:n_times){
    Cases<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]==cause)
    Controls_1<-(Mat_data[,"T"]> times[t] )
    Controls_2<-(Mat_data[,"T"]< times[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    if (weights$method!="marginal"){ 
      Weights_controls_1<-1/(weights$IPCW.times[,t]*n)  }
    else{
      Weights_controls_1<-rep(1/(weights$IPCW.times[t]*n),times=n)
    }
    Weights_controls_1<-Weights_controls_1[order_marker] 
    Weights_cases<-Weights_cases_all 
    Weights_controls_2<-Weights_cases_all
    Weights_cases[!Cases]<-0
    Weights_controls_1[!Controls_1]<-0
    Weights_controls_2[!Controls_2]<-0
    den_TP_t<-sum(Weights_cases)
    den_FP_1_t<-sum(Weights_controls_1)
    den_FP_2_t<-sum(Weights_controls_2)+sum(Weights_controls_1)
    if(den_TP_t!=0){  
      TP_tbis<-c(0,cumsum(Weights_cases))/den_TP_t
      TP_t<-TP_tbis[!duplicated(marker[order_marker])]
    }
    else TP_t<-NA
    if(den_FP_1_t!=0){
      FP_1_tbis<-c(0,cumsum(Weights_controls_1))/den_FP_1_t
      FP_1_t<-FP_1_tbis[!duplicated(marker[order_marker])]}
    else FP_1_t<-NA
    if(den_FP_2_t!=0){
      FP_2_tbis<-c(0,cumsum(Weights_controls_1)+cumsum(Weights_controls_2))/den_FP_2_t
      FP_2_t<-FP_2_tbis[!duplicated(marker[order_marker])]}
    else FP_2_t<-NA
    # internal fonction to compute an area under a curve by trapezoidal rule
    AireTrap<-function(Abs,Ord){
      nobs<-length(Abs)
      dAbs<-Abs[-1]-Abs[-nobs]
      mil<-(Ord[-nobs]+Ord[-1])/2
      area<-sum(dAbs*mil)
      return(area)
    }
    if ( den_TP_t*den_FP_1_t != 0){AUC_1[t]<-AireTrap(FP_1_t,TP_t)}
    else AUC_1[t]<-NA
    if ( den_TP_t*den_FP_2_t != 0){AUC_2[t]<-AireTrap(FP_2_t,TP_t)}
    else AUC_2[t]<-NA
    if(ROC==TRUE){ 
      TP[,t]<-TP_t
      FP_1[,t]<-FP_1_t
      FP_2[,t]<-FP_2_t
    }  
    CumInci[t]<-c(den_TP_t)
    surv[t]<-c(den_FP_1_t)
    Stats[t,]<-c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2))
  }
  # }}}
  inference<-NA
  if (iid==TRUE){   
    if(weighting!="marginal"){
      stop("Error : Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")
    } 
    else{      
      # create iid representation required for inference procedures
      out_iid<-vector("list", n_times)
      names(out_iid)<-paste("t=",times,sep="")
      vect_iid_comp_time<-rep(NA,times=n_times)
      names(vect_iid_comp_time)<-paste("t=",times,sep="")
      mat_iid_rep<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep)<-paste("t=",times,sep="")  
      mat_iid_rep_star<-matrix(NA,nrow=n,ncol=n_times)
      colnames(mat_iid_rep_star)<-paste("t=",times,sep="")  
      vetc_se<-rep(NA,times=n_times)
      names(vetc_se)<-paste("t=",times,sep="")
      vetc_sestar<-rep(NA,times=n_times)    
      names(vetc_sestar)<-paste("t=",times,sep="")   
     # compute iid for Kaplan Meier
     MatInt0TcidhatMCksurEff <- Compute.iid.KM(times=T,status=delta)
      for (j in 1:n_times){
        #compute iid representation when AUC can be computed
        if(!is.na(AUC_1[j]) | !is.na(AUC_2[j])){
          out_iid[[j]]<-compute_iid_decomposition(t=times[j],n=n,cause=cause,F01t=CumInci[j],St=surv[j],weights,T,delta,marker,MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)} 
        else{
          out_iid[[j]]<-NA}
        #browser()
        #save output for inference for AUC_1 when AUC_1 can be computed        
        if(!is.na(AUC_1[j])){
          mat_iid_rep_star[,j]<-out_iid[[j]]$iid_representation_AUCstar
          vetc_sestar[j]<-out_iid[[j]]$seAUCstar
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }
        #save output for inference for AUC_2 when AUC_2 can be computed         
        if(!is.na(AUC_2[j])){
          mat_iid_rep[,j]<-out_iid[[j]]$iid_representation_AUC
          vetc_se[j]<-out_iid[[j]]$seAUC
          vect_iid_comp_time[j]<-out_iid[[j]]$computation_times               
        }   
      }
      inference<-list(mat_iid_rep_2=mat_iid_rep,
                      mat_iid_rep_1=mat_iid_rep_star,
                      vect_sd_1=vetc_sestar,
                      vect_sd_2=vetc_se,
                      vect_iid_comp_time=vect_iid_comp_time
                      )
    }
  }
  stop_computation_time<-Sys.time() 
  # output if there is competing risks or not
  if (max(Stats[,3])==0){
    out <- list(TP=TP,FP=FP_1,AUC=AUC_1,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats[,c(1,2,4)],weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwsurvivalROC"
    out
  }else{
    out <- list(TP=TP,FP_1=FP_1,AUC_1=AUC_1,FP_2=FP_2,AUC_2=AUC_2,times=times,
                CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,weights=weights,
                inference=inference,computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),iid=iid)
    class(out) <- "ipcwcompetingrisksROC"
    out    
  }
}




