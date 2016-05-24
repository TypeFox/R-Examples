# {{{ Inputs :
# cutpoint       : the cutpoint for which we aim at estimating Se, Sp, PPV and NPV
# T              : vector of observed failure times
# delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
# marker         : vector ofmarker values
# other_markers  : (default is NULL, should be a matrix) other markers that can be associated with the censoring mechanism
# cause          : the value that indicates the main event of interest
# weighting      : (default is "marginal") weighting technique for IPCW : "marginal" for Kaplan-Meier, "cox" for proportional hazards cox model, "aalen" for additive aalen model
# times          : vector of times you want to compute the time dependent AUC.
# ROC            : if TRUE, then save True Positive fraction (Sensitivity) and False Positive fraction (1-Specificity)
#                  for all time in vetor times
#iid             : TRUE or FALSE, indicates if we want to compute the iid representation of the estimators
# }}}
# {{{ Outputs :
# TP
# FP_1
# FP_2
# PPV
# NPV_1
# NPV_2
# Stats              : matrix with numbers of observed cases and controls, and censored before time point, no. subjects with  maker > or < cutpoint
# inference          
# computation_time
# iid
# n
# times
# cutpoint
# }}}
SeSpPPVNPV <-function(cutpoint,T,delta,marker,other_markers=NULL,cause,weighting="marginal",times,iid=FALSE){
  ## browser()
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
  if(weighting!="marginal" & iid){
    stop("Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")
  } 
  # }}}
  # {{{ check if there are missing values, and delete rows with missing values
  if (weighting %in% c("cox","aalen") & !missing(other_markers) ){
    is_not_na <- as.logical(apply(!is.na(cbind(T,delta,marker,other_markers)),1,prod))
    T <- T[is_not_na]
    delta <- delta[is_not_na]
    marker <- marker[is_not_na]
    other_markers <- as.matrix(other_markers[is_not_na,])
  }else{
    is_not_na <- as.logical(apply(!is.na(cbind(T,delta,marker)),1,prod)) 
    T <- T[is_not_na]
    delta <- delta[is_not_na]
    marker <- marker[is_not_na]
  }
  # }}} 
  start_computation_time <- Sys.time()
  # {{{ create some usefull objects
  n <- length(T)
  n_times <- length(times)
  if (n_times==1){times<-c(0,times)
                  n_times<-2}           # trick to use ipcw.cox() even if there is only one time
  times <- times[order(times)]
  times_names <- paste("t=",times,sep="")
  # }}}
  # {{{ output initialisation
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  names(CumInci) <- times_names
  names(surv) <- times_names
  Stats <- matrix(NA,nrow=n_times,ncol=6)
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t","Positive (X>c)","Negative (X<=c)")
  rownames(Stats) <- times_names
  Stats[,c("Positive (X>c)","Negative (X<=c)")] <- matrix(c(sum(marker>cutpoint),sum(marker<=cutpoint)),byrow=TRUE,ncol=2,nrow=nrow(Stats))
  # }}}
  # {{{  computation of weights (1/2)
  # we need to order to use the ipcw() fonction
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
  ## order_marker<-order(-marker)
  ## Mat.data<-cbind(T,delta,marker)[order_marker,]
  Mat.data<-cbind(T,delta,marker)
  colnames(Mat.data)<-c("T","delta","marker")
  # Create some weights
  Weights_cases_all<-1/(weights$IPCW.subjectTimes*n)
  ## Weights_cases_all<-Weights_cases_all[order_marker]
  Weights_cases_all<-Weights_cases_all
  # }}}
  # {{{ Initialisation of outputs
  FP_1 <- rep(NA,n_times)
  TP <- rep(NA,n_times)
  FP_2 <- rep(NA,n_times)
  PPV <- rep(NA,n_times)
  NPV_1<- rep(NA,n_times)
  NPV_2<- rep(NA,n_times)
  names(FP_1) <- times_names
  names(TP) <- times_names
  names(FP_2) <- times_names
  names(PPV) <- times_names
  names(NPV_1) <- times_names
  names(NPV_2) <- times_names
  # }}}
  # {{{ if iid=TRUE, then compute KM censoring estimator iid representation sqrt(n)(G(T_i)-\hat{G}(T_i)), i.e at each observed time
  if (iid){
    MatInt0TcidhatMCksurEff <- Compute.iid.KM(T,delta!=0)
    epsiloni.Se <- matrix(NA,nrow=n,ncol=n_times)
    epsiloni.Sp1 <- matrix(NA,nrow=n,ncol=n_times)
    epsiloni.Sp2 <- matrix(NA,nrow=n,ncol=n_times)
    epsiloni.PPV <- matrix(NA,nrow=n,ncol=n_times)
    epsiloni.NPV1 <- matrix(NA,nrow=n,ncol=n_times)
    epsiloni.NPV2 <- matrix(NA,nrow=n,ncol=n_times)   
    se.Se <- rep(NA,n_times)
    se.Sp1 <- rep(NA,n_times)
    se.Sp2 <- rep(NA,n_times)
    names(se.Sp1) <- times_names
    names(se.Sp2) <- times_names
    names(se.Se) <- times_names
    colnames(epsiloni.Se) <- times_names
    colnames(epsiloni.Sp1) <- times_names
    colnames(epsiloni.Sp2) <- times_names
    colnames(epsiloni.PPV) <- times_names
    colnames(epsiloni.NPV1) <- times_names
    colnames(epsiloni.NPV2) <- times_names
  }
  # }}}
  ## browser()
  for(t in 1:n_times){
    # browser()
    # {{{ first part of loop on all timepoints t
    Cases<-(Mat.data[,"T"]< times[t] &  Mat.data[,"delta"]==cause)
    Controls_1<-(Mat.data[,"T"]> times[t] )
    Controls_2<-(Mat.data[,"T"]< times[t] &  Mat.data[,"delta"]!=cause & Mat.data[,"delta"]!=0)  
    if (weights$method!="marginal"){ 
      Weights_controls_1<-1/(weights$IPCW.times[,t]*n)
    }else{
      Weights_controls_1<-rep(1/(weights$IPCW.times[t]*n),times=n)
    }
    ## Weights_controls_1<-Weights_controls_1[order_marker]
    Weights_controls_1<-Weights_controls_1
    Weights_cases<-Weights_cases_all 
    Weights_controls_2<-Weights_cases_all
    Weights_cases[!Cases]<-0
    Weights_controls_1[!Controls_1]<-0
    Weights_controls_2[!Controls_2]<-0
    den_TP_t<-sum(Weights_cases)
    den_FP_1_t<-sum(Weights_controls_1)
    den_FP_2_t<-sum(Weights_controls_2)+sum(Weights_controls_1)
    # compute P(X>c) and P(X<c) with weighted sums
    den_PPV_t <- sum(Weights_cases[which(Mat.data[,"marker"] > cutpoint)] + Weights_controls_1[which(Mat.data[,"marker"] > cutpoint)] + Weights_controls_2[which(Mat.data[,"marker"] > cutpoint)]) 
    den_NPV_t <- sum(Weights_cases[which(Mat.data[,"marker"] <= cutpoint)] + Weights_controls_1[which(Mat.data[,"marker"] <= cutpoint)] + Weights_controls_2[which(Mat.data[,"marker"] <= cutpoint)]) 
    if(den_TP_t!=0){  
      TP[t] <- sum(Weights_cases[which(Mat.data[,"marker"] > cutpoint)])/den_TP_t
    }
    if(den_FP_1_t!=0){
      FP_1[t] <- sum(Weights_controls_1[which(Mat.data[,"marker"] > cutpoint)])/den_FP_1_t
    }
    if(den_FP_2_t!=0){
      FP_2[t] <- sum(Weights_controls_1[which(Mat.data[,"marker"] > cutpoint)] + Weights_controls_2[which(Mat.data[,"marker"] > cutpoint)])/den_FP_2_t
      NPV_2[t] <- ((1-FP_2[t])*den_FP_2_t)/den_NPV_t
    }
    # compute/save decriptive statistics
    CumInci[t] <- c(den_TP_t)
    surv[t] <- c(den_FP_1_t)
    Stats[t,1:4] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2))
    # compute PPV and NPV 1
    PPV[t] <- (TP[t]*CumInci[t])/den_PPV_t
    NPV_1[t] <- ((1-FP_1[t])*surv[t])/den_NPV_t
    # }}}
    # {{{ iid representtion
    if (iid){
      ## browser()
      Int0tdMCsurEffARisk <- NA
      # compute iid sqrt(n)(G(t)-\hat{G}(t)), i.e at time point
      Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[max(which(Mat.data[,"T"] <= times[t])),]
      Weights_cases <- Weights_cases
      Weights_controls_2 <- Weights_controls_2
      Weights_controls_1 <- Weights_controls_1
      epsiloni.Se[,t] <- Weights_cases*n*(as.numeric(Mat.data[,"marker"] > cutpoint)-TP[t])/CumInci[t] + colMeans(MatInt0TcidhatMCksurEff*(Weights_cases*n*(as.numeric(Mat.data[,"marker"] > cutpoint)-TP[t])/CumInci[t]))   
      epsiloni.Sp1[,t] <- (n/sum(Controls_1))*as.numeric(Mat.data[,"T"]>times[t])*(as.numeric(Mat.data[,"marker"] <= cutpoint)-(1-FP_1[t]))
      epsiloni.Sp2[,t] <- (Weights_controls_2*n + Weights_controls_1*n)*(as.numeric(Mat.data[,"marker"] <= cutpoint)-(1-FP_2[t]))/(1-CumInci[t]) + colMeans(MatInt0TcidhatMCksurEff*((Weights_controls_2*n)*(as.numeric(Mat.data[,"marker"] <= cutpoint)-(1-FP_2[t]))/(1-CumInci[t]))) + mean((Weights_controls_1*n)*(as.numeric(Mat.data[,"marker"] <= cutpoint)-(1-FP_2[t]))/(1-CumInci[t]))*Int0tdMCsurEffARisk    
      # To check iid representation of Sp1
      # sd(epsiloni.Sp1)/sqrt(n)*sqrt((n-1)/n)
      # sqrt(Sp1*(1-Sp1)/sum(Controls.1))
      epsiloni.PPV[,t] <- (as.numeric(Mat.data[,"marker"] > cutpoint)/den_PPV_t)*((Weights_cases*n + Weights_controls_2*n)*(as.numeric(Mat.data[,"delta"]==cause)-as.numeric(Mat.data[,"delta"]!=0)*PPV[t]
                                                                                                                            ) - Weights_controls_1*n*PPV[t]
                                                                                  ) + colMeans(MatInt0TcidhatMCksurEff*((as.numeric(Mat.data[,"marker"] > cutpoint)/den_PPV_t
                                                                                                                         )*(Weights_cases*n + Weights_controls_2*n
                                                                                                                            )*(as.numeric(Mat.data[,"delta"]==cause)-as.numeric(Mat.data[,"delta"]!=0)*PPV[t]
                                                                                                                               ))
                                                                                               )-mean((as.numeric(Mat.data[,"marker"] > cutpoint)/den_PPV_t
                                                                                                       )* Weights_controls_1*n*PPV[t]
                                                                                                      )*Int0tdMCsurEffARisk    
      epsiloni.NPV2[,t] <- (as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t)*((Weights_cases*n + Weights_controls_2*n)*(as.numeric(Mat.data[,"delta"]!=0 & Mat.data[,"delta"]!=cause)-as.numeric(Mat.data[,"delta"]!=0)*NPV_2[t]
                                                 ) + Weights_controls_1*n*(1-NPV_2[t])
                                                 ) + colMeans(MatInt0TcidhatMCksurEff*((as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t
                                                                                        )*(Weights_cases*n + Weights_controls_2*n
                                                                                           )*(as.numeric(Mat.data[,"delta"]!=0 & Mat.data[,"delta"]!=cause)-as.numeric(Mat.data[,"delta"]!=0)*NPV_2[t]
                                                                                              ))
                                                              )+mean((as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t
                                                                      )* Weights_controls_1*n*(1-NPV_2[t]))*Int0tdMCsurEffARisk    
      epsiloni.NPV1[,t] <- (as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t)*((Weights_cases*n + Weights_controls_2*n)*(-as.numeric(Mat.data[,"delta"]!=0)*NPV_1[t]
                                                 ) + Weights_controls_1*n*(1-NPV_1[t])
                                                 )  + colMeans(MatInt0TcidhatMCksurEff* ((as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t
                                                                                          )*(Weights_cases*n + Weights_controls_2*n
                                                                                             )*(-as.numeric(Mat.data[,"delta"]!=0)*NPV_1[t]
                                                                                                ))
                                                               )+mean((as.numeric(Mat.data[,"marker"] <= cutpoint)/den_NPV_t
                                                                       )*Weights_controls_1*n*(1-NPV_1[t])
                                                                      )*Int0tdMCsurEffARisk
    }
    # }}}
  }
  ## browser()
  # {{{ compute se
  if (iid){
    se.Se <- apply(epsiloni.Se,2,sd)/sqrt(n)
    se.Sp1 <- apply(epsiloni.Sp1,2,sd)/sqrt(n)
    se.Sp2 <- apply(epsiloni.Sp2,2,sd)/sqrt(n)
    se.PPV <- apply(epsiloni.PPV,2,sd)/sqrt(n)
    se.NPV1 <- apply(epsiloni.NPV1,2,sd)/sqrt(n)
    se.NPV2 <- apply(epsiloni.NPV2,2,sd)/sqrt(n)
  }
  ## browser()
  # }}}
  # {{{ if iid then  build inference output
  if (iid){
    inference<-list(mat_iid_rep_Se=epsiloni.Se,
                    mat_iid_rep_Sp1=epsiloni.Sp1,
                    mat_iid_rep_Sp2=epsiloni.Sp2,
                    mat_iid_rep_PPV=epsiloni.PPV,
                    mat_iid_rep_NPV1=epsiloni.NPV1,
                    mat_iid_rep_NPV2=epsiloni.NPV2,
                    vect_se_Se=se.Se,
                    vect_se_Sp1=se.Sp1,
                    vect_se_Sp2=se.Sp2,
                    vect_se_PPV=se.PPV,
                    vect_se_NPV1=se.NPV1,
                    vect_se_NPV2=se.NPV2
                    )
  }else{
    inference<-NA
  }
  # }}}
  ## browser()
  stop_computation_time<-Sys.time()
  # {{{ output if there is competing risks or not
  if (max(Stats[,3])==0){
    out <- list(TP=TP,FP=FP_1,PPV=PPV,NPV=NPV_1,Stats=Stats[,-3],
                inference=inference,
                computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),
                iid=iid,
                n=n,times=times,weights=weights,
                cutpoint=cutpoint)
    class(out) <- "ipcwsurvivalSeSpPPVNPV"
    out
  }else{
    out <- list(TP=TP,FP_1=FP_1,FP_2=FP_2,PPV=PPV,NPV_1=NPV_1,NPV_2=NPV_2,
                Stats=Stats,
                inference=inference,
                computation_time=difftime(stop_computation_time,start_computation_time,units="secs"),
                iid=iid,
                n=n,times=times,weights=weights,
                cutpoint=cutpoint)
    class(out) <- "ipcwcompetingrisksSeSpPPVNPV"
    out    
  }
  # }}}
}
