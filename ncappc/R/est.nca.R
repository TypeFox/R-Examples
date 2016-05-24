# R function to estimate NCA metrics
# Chayan, 12/2014

# roxygen comments
#' Estimates individual NCA metrics from concentration vs time data.
#' data.
#'
#' \pkg{est.nca} estimates a comprehensive set of NCA metrics
#' for a given individual using concentration vs time profile.
#'
#' \pkg{est.nca} estimates a comprehensive set of NCA metrics using the
#' concentration-time profile of an individual. NCA metrics are eatimated
#' according to traditional PK calculations. The names of the various NCA
#' metrics estimated in this package are assigned mainly following the names
#' used in WinNonlin. This package accepts any of the three different types of
#' drug administration, (i) iv-bolus, (ii) iv-infusion and (iii) extravascular;
#' \pkg{ncappc} also can accept both non-steady state and steady-state data. 
#' The NCA metrics that are estimated and reported by \pkg{ncappc} are listed 
#' below.
#' \itemize{
#'  \item \strong{C0} is the initial concentration at the dosing time. It is 
#'  the observed concentration at the dosing time, if available. Otherwise 
#'  it is approximated using the following rules.
#'  \item \strong{Cmax, Tmax and Cmax_D} are the value and the time of maximum
#'  observed concentration, respectively. If the maximum concentration is not
#'  unique, the first maximum is used. For steady state data, The maximum value
#'  between the dosing intervals is considered. Cmax_D is the dose normalized
#'  maximum observed concentration.
#'  \item \strong{Clast and Tlast} are the last measurable positive
#'  comcentration and the corresponding time, respectively. 
#'  \item \strong{AUClast} is the area under the concentration vs. time curve
#'  from the first observed to last measurable concentration.
#'  \item \strong{AUMClast} is the first moment of the concentration vs. time
#'  curve from the first observed to last measurable concentration.
#'  \item \strong{MRTlast} is the mean residence time from the first observed 
#'  to last measurable concentration.
#'  \item \strong{No_points_Lambda_z} is the number of observed data points 
#'  used to determine the best fitting regression line in the elimination 
#'  phase.
#'  \item \strong{AUC_pBack_Ext_obs} is the percentage of AUCINF_obs that is 
#'  contributed by the back extrapolation to estimate C0.
#'  \item \strong{AUC_pBack_Ext_pred} is the percentage of AUCINF_pred that is 
#'  contributed by the back extrapolation to estimate C0.
#'  \item \strong{AUClower_upper} is the AUC under the concentration-time 
#'  profile within the user-specified window of time privided as the 
#'  "AUCTimeRange" argument. In case of empty "AUCTimeRange" argument, 
#'  AUClower_upper is equal to the AUClast.
#'  \item \strong{Rsq, Rsq_adjusted and Corr_XY} are regression coefficient 
#'  of the regression line used to estimate the elimination rate constant, the 
#'  adjusted value of Rsq and the square root of Rsq, respectively.
#'  \item \strong{Lambda_z} is the elimination rate constant estimated from the 
#'  regression line representing the terminal phase of the concentration-time 
#'  prifile.
#'  \item \strong{Lambda_lower and Lambda_upper} are the lower and upper limit 
#'  of the time values from the concentration-time profile used to estimate 
#'  Lambda_z, respectively, in case the "LambdaTimeRange" is used to specify 
#'  the time range.
#'  \item \strong{HL_Lambda_z} is terminal half-life of the drug.
#'  \item \strong{AUCINF_obs and AUCINF_obs_D} are AUC estimated from the first 
#'  sampled data extrapolated to infinity and its dose normalized version, 
#'  respectively. The extrapolation in the terminal phase is based on the last 
#'  observed concentration Clast_obs.
#'  \item \strong{AUC_pExtrap_obs} is the percentage of the AUCINF_obs that is 
#'  contributed by the extrapolation from the last sampling time to infinity.
#'  \item \strong{AUMCINF_obs} is AUMC estimated from the first sampled data 
#'  extrapolated to infinity. The extrapolation in the terminal phase is based 
#'  on the last observed concentration.
#'  \item \strong{AUMC_pExtrap_obs} is the percentage of the AUMCINF_obs that 
#'  is contributed by the extrapolation from the last sampling time to infinity.
#'  \item \strong{Vz_obs} is the volume of distribution estimated based on 
#'  total AUC
#'  \item \strong{Cl_obs} is total body clearance.
#'  \item \strong{AUCINF_pred and AUCINF_pred_D} are AUC from the first sampled 
#'  data extrapolated to infinity and its dose normalized version, respectively.
#'  The extrapolation in the terminal phase is based on the last predicted 
#'  concentration obtained from the regression line used to estimate Lambda_z 
#'  (Clast_pred).
#'  \item \strong{AUC_pExtrap_pred} is the percentage of the AUCINF_pred that 
#'  is contributed by the extrapolation from the last sampling time to infinity.
#'  \item \strong{AUMCINF_pred} is AUMC estimated from the first sampled data 
#'  extrapolated to infinity. The extrapolation in the terminal phase is based 
#'  on the last predicted concentration obtained from the regression line used 
#'  to estimate Lambda_z (Clast_pred).
#'  \item \strong{AUMC_pExtrap_pred} is the percentage of the AUMCINF_pred that 
#'  is contributed by the extrapolation from the last sampling time to infinity.
#'  \item \strong{Vz_pred} is the volume of distribution estimated based on 
#'  AUCINF_pred.
#'  \item \strong{Cl_pred} is the total body clearance estimated based on 
#'  AUCINF_pred.
#'  \item \strong{MRTINF_obs} is the mean residence time from the first 
#'  sampled time extrapolated to infinity based on the last observed 
#'  concentration (Clast_obs).
#'  \item \strong{MRTINF_pred} is the mean residence time from the first 
#'  sampled time extrapolated to infinity based on the last predicted 
#'  concentration obtained from the regression line used to estimate Lambda_z 
#'  (Clast_pred).
#'  \item \strong{Tau} is the dosing interval for steady-state data. This value 
#'  is assumed equarion over multiple doses.
#'  \item \strong{Cmin and Tmin} are the minimum concentration between 0 and 
#'  Tau and the corresponding time, respectively.
#'  \item \strong{Cavg} is the average concentration between 0 and Tau for 
#'  steady-state data.
#'  \item \strong{AUCtau and AUMCtau} are AUC and AUMC between 0 and Tau for 
#'  steady-state data.
#'  \item \strong{Clss} is an estimate of the total body clearance for
#'  steady-state data.
#'  \item \strong{Vss_obs and Vss_pred} are estimated volume of distribution at 
#'  steady-state based on Clast_obs and Clast_pred, respectively.
#'  \item \strong{p_Fluctuation} is the percentage of the fluctuation of the 
#'  concentration between 0 and Tau for steady-state data.
#'  \item \strong{Accumulation_Index} is \eqn{1/(1-e^(-\lambda_z*\tau))}
#' }
#' 
#' @param time Numeric array for time
#' @param conc Numeric array for concentration
#' @param backExtrp If back-extrapolation is needed for AUC (TRUE or FALSE)
#'   (\strong{FALSE})
#' @param negConcExcl Exclude -ve conc (\strong{FALSE})
#' @param doseType Steady-state (ss) or nonsteady-state (ns) dose
#'   (\strong{"ns"})
#' @param adminType Route of administration
#'   (iv-bolus,iv-infusion,extravascular) (\strong{"extravascular"})
#' @param doseAmt Dose amounts (\strong{"NULL"})
#' @param method linear, log or linear-log (\strong{"linear-log"})
#' @param AUCTimeRange User-defined window of time used to estimate AUC 
#'   (\strong{"NULL"})
#' @param LambdaTimeRange User-defined window of time to estimate elimination 
#'   rate-constant (\strong{"NULL"})
#' @param LambdaExclude User-defined excluded observation time points for
#'   estimation of elimination rate-constant (\strong{"NULL"})
#' @param Tau Dosing interval for steady-state data (\strong{"NULL"})
#' @param doseTime Dose time prior to the first observation for steady-state data (\strong{\code{NULL}})
#' @param TI Infusion duration (\strong{"NULL"})
#' @param simFile Name of the simulated concentration-time data if present
#'   (\strong{"NULL"})
#' @param dset Type, i.e., observed or simulated concentration-time data set
#'   ("obs" or "sim") (\strong{"obs"})
#' 
#' @return An array of estimated NCA metrics
#' @export
#'


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# function to estimate NCA parameters

est.nca <- function(time,conc,backExtrp=FALSE,negConcExcl=FALSE,doseType="ns",adminType="extravascular",
                    doseAmt=NULL,method="linear-log",AUCTimeRange=NULL,LambdaTimeRange=NULL,LambdaExclude=NULL,
                    doseTime=doseTime,Tau=NULL,TI=NULL,simFile=NULL,dset="obs"){

  "tail" <- "head" <- "lm" <- "coef" <- NULL
  rm(list=c("tail","head","lm","coef"))
  
  ## set all NCA metrics, except Tau to "NaN"
  #C0 <- "NaN"; Tmax <- "NaN"; Cmax <- "NaN"; Cmax_D <- "NaN"; Tlast <- "NaN"; Clast <- "NaN"; AUClast <- "NaN"; AUMClast <- "NaN"; MRTlast <- "NaN"; No_points_Lambda_z <- "NaN"
  #AUC_pBack_Ext_obs <- AUC_pBack_Ext_pred <- "NaN"; AUClower_upper <- "NaN"; Rsq <- "NaN"; Rsq_adjusted <- "NaN"; Corr_XY <- "NaN"; Lambda_z <- "NaN"; Lambda_z_lower <- "NaN"; Lambda_z_upper <- "NaN"; HL_Lambda_z <- "NaN"
  #AUCINF_obs <- "NaN"; AUCINF_D_obs <- "NaN"; AUC_pExtrap_obs <- "NaN"; Vz_obs <- "NaN"; Cl_obs <- "NaN"; AUCINF_pred <- "NaN"; AUCINF_D_pred <- "NaN"; AUC_pExtrap_pred <- "NaN"; Vz_pred <- "NaN"; Cl_pred <- "NaN"
  #AUMCINF_obs <- "NaN"; AUMC_pExtrap_obs <- "NaN"; AUMCINF_pred <- "NaN"; AUMC_pExtrap_pred <- "NaN"; MRTINF_obs <- "NaN"; MRTINF_pred <- "NaN"; Tmin <- "NaN"; Cmin <- "NaN"; Cavg <- "NaN"
  #AUCtau <- "NaN"; AUMCtau <- "NaN"; Clss <- "NaN"; Vss_obs <- "NaN"; Vss_pred <- "NaN"; p_Fluctuation <- "NaN"; Accumulation_Index <- "NaN"
  
  Tau <- ifelse(doseType=="ns", NA, Tau)
  
  ntime <- as.numeric(time)
  nconc <- as.numeric(conc)
  
  # Exclude zero or negative concentration from calculations
  if (negConcExcl == TRUE){
    zid <- which (nconc < 0)
    if (length(zid) > 0){ntime  <- ntime[-zid]; nconc  <- nconc[-zid]}
  }
  
  if (length(nconc) >= 2){
    # Steady state data
    if (doseType == "ss"){
      if(is.null(doseTime)){
        # Index for the first observation with non-zero conc
        ssIdx1 <- min(which(nconc>0))
        # Steady state 1st observation time and conc
        sc1 <- nconc[ssIdx1]; st1 <- ntime[ssIdx1]
        # Steady state last observation time and conc
        ssIdx2 <- which.min(abs(ntime-(st1+Tau)))  # Index for the last ss time (Tau added to first time and find the nearest time to that number)
        ssIdx2 <- max(which(abs(round(ntime,0)-lTm)==0))  # Find the index in time vector that matches
        sc2 <- nconc[ssIdx2]; st2 <- ntime[ssIdx2]
      }else{
        # Index for the first observation within the steady state interval
        ssIdx1 <- min(which(abs(round(ntime,0)-doseTime)==0))
        # Steady state 1st observation time and conc
        sc1 <- nconc[ssIdx1]; st1 <- ntime[ssIdx1]
        # Steady state last observation time and conc
        lTm <- round(ntime[ssIdx1],0)+Tau  # Tau added to first time
        ssIdx2 <- which(ntime==(st1+Tau-min(abs(round(ntime,0)-lTm))[1]))  # Find the index in time vector closest to first obs time + Tau
        #ssIdx2 <- max(which(abs(round(ntime,0)-lTm)==0))  # Find the index in time vector that matches
        sc2 <- nconc[ssIdx2]; st2 <- ntime[ssIdx2]
      }
      # steady state observations
      nconc <- nconc[which(ntime>=st1 & ntime<=st2)]
      ntime <- ntime[which(ntime>=st1 & ntime<=st2)]
      otime <- ntime - min(ntime)       # Time if off-set to zero for estimation of the 1st moment
      ssnPt <- length(ntime)
      
      #sconc <- nconc[which(ntime>=0 & ntime<=Tau)]; stime <- ntime[which(ntime>=0 & ntime<=Tau)]; ssnPt <- length(stime)
      if (ssnPt != 0){
        mIdx <- which(nconc == min(nconc[nconc!=0]))[1]; Tmin <- ntime[mIdx]; Cmin <- nconc[mIdx]
        AUCtau <- AUMCtau <- 0.
      }
    }
    
    # Calculation of C0
    noBackTime <- TRUE
    if (backExtrp == TRUE){
      if (ntime[1] == 0){
        C0 <- as.numeric(nconc[1])
      }else{
        if (adminType == "extravascular" | adminType == "iv-infusion"){
          C0 <- 0
          nconc <- c(C0, nconc)
          if (doseType == "ss" && ssnPt != 0) if (ssnPt != 0){C0 <- as.numeric(min(nconc))}
        }else if (adminType == "iv-bolus"){
          slope <- (nconc[2]-nconc[1])/(ntime[2]-ntime[1])
          if (nconc[1]==0 | nconc[2]==0){
            C0 <- as.numeric(nconc[nconc>0][1])
          }else if (slope >= 0){
            C0 <- as.numeric(nconc[1])
          }else{
            C0 <- as.numeric(exp((ntime[1]*log(nconc[2])-ntime[2]*log(nconc[1]))/(ntime[1]-ntime[2])))
            nconc <- c(C0, nconc) # extrapolation via log-linear regression
            ntime <- c(0, ntime)
            noBackTime <- FALSE
          }
        }
      }
    }
    
    # AUClast & AUMClast and parameters for steady state data
    otime   <- ntime - min(ntime)   # Time if off-set to zero for estimation of the 1st moment
    AUClast <- AUMClast <- AUCnoC0 <- AUClower_upper <- 0.
    
    nPt    <- length(nconc)   # No. of data points
    Cmax   <- max(nconc)
    mxId   <- which(nconc == Cmax)[length(which(nconc == Cmax))]
    Tmax   <- ntime[mxId]
    Cmax_D <- ifelse(!is.null(doseAmt), nconc[mxId]/doseAmt, NA)
    lIdx   <- max(which(nconc == tail(nconc[nconc>0],1))) # Index for last positive concentration
    Tlast  <- ntime[lIdx]; Clast <- nconc[lIdx]
    
    for(r in 1:(nPt-1)){
      if (method == "linear"){
        delauc   <- mean(nconc[r:(r+1)])*(ntime[r+1]-ntime[r])
        delaumc  <- 0.5*((nconc[r+1]*otime[r+1])+(nconc[r]*otime[r]))*(otime[r+1]-otime[r])
        AUClast  <- sum(AUClast, delauc)
        AUMClast <- sum(AUMClast, delaumc)
        if (backExtrp == TRUE){
          if (noBackTime == TRUE){
            AUCnoC0 <- sum(AUCnoC0, delauc)
          }else if (noBackTime == FALSE & r>1){
            AUCnoC0 <- sum(AUCnoC0, delauc)
          }
        }
        if (!is.null(AUCTimeRange)){if (ntime[r] >= min(AUCTimeRange) & ntime[r+1] <= max(AUCTimeRange)){AUClower_upper <- sum(AUClower_upper, delauc)}}
        if (doseType == "ss" && ssnPt != 0){
          AUCtau  <- sum(AUCtau, delauc); AUMCtau <- sum(AUMCtau, delaumc)
        }
      }else if (method == "log"){
        delauc   <- (nconc[r+1]-nconc[r])*(ntime[r+1]-ntime[r])/log(nconc[r+1]/nconc[r])
        delaumc  <- ((((nconc[r+1]*otime[r+1])-(nconc[r]*otime[r]))*(otime[r+1]-otime[r]))/(log(nconc[r+1]/nconc[r]))) - ((nconc[r+1]-nconc[r])*((otime[r+1]-otime[r])**2)/((log(nconc[r+1]/nconc[r]))**2))
        AUClast  <- sum(AUClast, delauc)
        AUMClast <- sum(AUMClast, delaumc)
        if (backExtrp == TRUE){
          if (noBackTime == TRUE){
            AUCnoC0 <- sum(AUCnoC0, delauc)
          }else if (noBackTime == FALSE & r>1){
            AUCnoC0 <- sum(AUCnoC0, delauc)
          }
        }
        if (!is.null(AUCTimeRange)){if (ntime[r] >= min(AUCTimeRange) & ntime[r+1] <= max(AUCTimeRange)){AUClower_upper <- sum(AUClower_upper, delauc)}}
        if (doseType == "ss" && ssnPt != 0){
          AUCtau <- sum(AUCtau, delauc); AUMCtau <- sum(AUMCtau, delaumc)
        }
      }else if (method == "linear-log"){
        if (nconc[r+1]>=nconc[r] | nconc[r+1]<=0 | nconc[r]<=0){
          delauc   <- mean(nconc[r:(r+1)])*(ntime[r+1]-ntime[r])
          delaumc  <- 0.5*((nconc[r+1]*otime[r+1])+(nconc[r]*otime[r]))*(otime[r+1]-otime[r])
          AUClast  <- sum(AUClast, delauc)
          AUMClast <- sum(AUMClast, delaumc)
          if (backExtrp == TRUE){
            if (noBackTime == TRUE){
              AUCnoC0 <- sum(AUCnoC0, delauc)
            }else if (noBackTime == FALSE & r>1){
              AUCnoC0 <- sum(AUCnoC0, delauc)
            }
          }
          if (!is.null(AUCTimeRange)){if (ntime[r] >= min(AUCTimeRange) & ntime[r+1] <= max(AUCTimeRange)){AUClower_upper <- sum(AUClower_upper, delauc)}}
          if (doseType == "ss" && ssnPt != 0){
            AUCtau <- sum(AUCtau, delauc); AUMCtau <- sum(AUMCtau, delaumc)
          }
        }else{
          delauc   <- (nconc[r+1]-nconc[r])*(ntime[r+1]-ntime[r])/log(nconc[r+1]/nconc[r])
          delaumc  <- ((((nconc[r+1]*otime[r+1])-(nconc[r]*otime[r]))*(otime[r+1]-otime[r]))/(log(nconc[r+1]/nconc[r]))) - ((nconc[r+1]-nconc[r])*((otime[r+1]-otime[r])**2)/((log(nconc[r+1]/nconc[r]))**2))
          AUClast  <- sum(AUClast, delauc)
          AUMClast <- sum(AUMClast, delaumc)
          if (backExtrp == TRUE){
            if (noBackTime == TRUE){
              AUCnoC0 <- sum(AUCnoC0, delauc)
            }else if (noBackTime == FALSE & r>1){
              AUCnoC0 <- sum(AUCnoC0, delauc)
            }
          }
          if (!is.null(AUCTimeRange)){if (ntime[r] >= min(AUCTimeRange) & ntime[r+1] <= max(AUCTimeRange)){AUClower_upper <- sum(AUClower_upper, delauc)}}
          if (doseType == "ss" && ssnPt != 0){
            AUCtau <- sum(AUCtau, delauc); AUMCtau <- sum(AUMCtau, delaumc)
          }
        }
      }
    }
    if (is.null(AUCTimeRange)){AUClower_upper <- AUClast}
    if (backExtrp == TRUE){
      AUC_pBack_Ext_obs <- 100*(AUClast-AUCnoC0)/AUClast
      AUC_pBack_Ext_pred <- AUC_pBack_Ext_obs
    }
    
    # MRTlast
    if (adminType != "iv-infusion"){
      MRTlast <- ifelse(doseType != "ss", AUMClast/AUClast, AUMCtau/AUCtau)
    }else{
      MRTlast <- (AUMClast/AUClast)-(TI/2)
    }
    
    # NCA metrics involving terminal-phase elimination
    # Determine the lower and upper indeces of the time vector used for Lambda calculation
    AUCINF_obs <- AUCINF_pred <- AUMCINF_obs <- AUMCINF_pred <- 0.
    
    llower <- tail(which(nconc == max(nconc)),1); lupper <- nPt
    
    # Exclude Cmax from elimination phase extrapolation for non-bolus dose
    if(adminType!="iv-bolus") llower <- llower+1
    
    if (!is.null(LambdaTimeRange)){
      Lambda_z_lower <- sort(LambdaTimeRange)[1]; Lambda_z_upper <- sort(LambdaTimeRange)[2]
      if (length(ntime[ntime >= Lambda_z_lower]) != 0 & length(ntime[ntime <= Lambda_z_upper]) != 0){
        llower <- head(match(ntime[ntime >= Lambda_z_lower]), 1); lupper <- tail(match(ntime[ntime <= Lambda_z_upper]), 1)
      }
    } #else{Lambda_z_lower <- "NaN"; Lambda_z_upper <- "NaN"}
    
    lconc <- nconc[llower:lupper]; ltime <- ntime[llower:lupper]
    # exclude time points as mentioned by the user
    if (!is.null(LambdaExclude)){lconc <- lconc[!ltime%in%LambdaExclude]; ltime <- ltime[!ltime%in%LambdaExclude]}
    # exclude points with zero or negative concentration
    zid <- which(lconc <= 0); if (length(zid) > 0){lconc <- lconc[-zid]; ltime <- ltime[-zid]}
    lconc <- log(lconc)
    lnPt  <- length(lconc)        # No. of data points for the descending part
    if (lnPt >= 3){
      # Determine the log-linear regression coefficients to calculate Lambda
      infd <- data.frame(np=numeric(0),rsq=numeric(0),arsq=numeric(0),m=numeric(0),inpt=numeric(0))
      for (r in 1:(lnPt-2)){
        mlr  <- lm(lconc[r:lnPt]~ltime[r:lnPt])
        if (is.na(coef(mlr)[2])) next
        if (coef(mlr)[2] >= 0) next
        n    <- (lnPt-r)+1
        rsq  <- summary(mlr)$r.squared
        trsq <- 1-((1-rsq)*(n-1)/(n-2))  # adjusted r^2
        tmp  <- cbind(np=n,rsq=rsq,arsq=trsq,m=(coef(mlr)[2]),inpt=(coef(mlr)[1]))
        infd <- rbind(infd,tmp)
      }
      if (nrow(infd) != 0){
        infd <- infd[order(infd$arsq,decreasing=T),]
        Rsq  <- infd$rsq[1]; Rsq_adjusted <- infd$arsq[1]; No_points_Lambda_z <- infd$np[1]; slope <- infd$m[1]; intercept <- infd$inpt[1]
        if (nrow(infd) > 1){
          for (r in 2:nrow(infd)){
            tarsq <- infd$arsq[r]; tDPt <- infd$np[r]
            if (Rsq_adjusted - tarsq > 0.0001) break
            if (No_points_Lambda_z > tDPt) next
            Rsq <- infd$rsq[r]; Rsq_adjusted <- tarsq; No_points_Lambda_z <- tDPt; slope <- infd$m[r]; intercept <- infd$inpt[r]
          }
        }
        Corr_XY      <- -1*sqrt(Rsq)
        Lambda_z     <- (-1*slope)
        HL_Lambda_z  <- log(2)/Lambda_z
        AUCINF_obs   <- exp(lconc[lnPt])/Lambda_z
        AUMCINF_obs  <- (exp(lconc[lnPt])/(Lambda_z**2))+(ltime[lnPt]*exp(lconc[lnPt])/Lambda_z)
        lastPt       <- exp((slope*ltime[lnPt])+intercept)
        AUCINF_pred  <- lastPt/Lambda_z
        AUMCINF_pred <- (lastPt/(Lambda_z**2))+(ltime[lnPt]*lastPt/Lambda_z)
        
        if(!is.null(AUCTimeRange) && (max(AUCTimeRange)>max(ltime))){
          uptime  <- max(AUCTimeRange)
          lowtime <- ifelse(min(AUCTimeRange) < max(ltime), max(ltime), min(AUCTimeRange))
          upconc  <- ifelse(exp(slope*uptime+intercept)==0, 0.0001, exp(slope*uptime+intercept))
          lowconc <- ifelse(lowtime==max(ltime), exp(lconc[lnPt]), exp(slope*lowtime+intercept))
          delauc  <- (upconc-lowconc)*(uptime-lowtime)/log(upconc/lowconc)
          AUClower_upper <- AUClower_upper + delauc
        }
        
        if (doseType == "ss" && (Tau > max(ltime) & AUCtau != 0)){
          uptime  <- Tau
          lowtime <- max(ltime)
          upconc  <- ifelse(exp(slope*uptime+intercept)==0, 0.0001, exp(slope*uptime+intercept))
          lowconc <- exp(lconc[lnPt])
          delauc  <- (upconc-lowconc)*(uptime-lowtime)/log(upconc/lowconc)
          delaumc <- ((((upconc*uptime)-(lowconc*lowtime))*(uptime-lowtime))/(log(upconc/lowconc))) - ((upconc-lowconc)*((uptime-lowtime)**2)/((log(upconc/lowconc))**2))
          AUCtau  <- sum(AUCtau, delauc)
          AUMCtau <- sum(AUMCtau, delaumc)
        }
      }
      
      if (AUClast != 0 & AUCINF_obs != 0){
        AUCINF_obs  <- AUClast+AUCINF_obs;  AUCINF_D_obs  <- ifelse(!is.null(doseAmt),AUCINF_obs/doseAmt,NA);  AUC_pExtrap_obs  <- 100*(AUCINF_obs-AUClast)/AUCINF_obs
        AUCINF_pred <- AUClast+AUCINF_pred; AUCINF_D_pred <- ifelse(!is.null(doseAmt),AUCINF_pred/doseAmt,NA); AUC_pExtrap_pred <- 100*(AUCINF_pred-AUClast)/AUCINF_pred
        if(doseType=="ns"){
          Vz_obs  <- ifelse(!is.null(doseAmt),doseAmt/(Lambda_z*AUCINF_obs),NA);  Cl_obs  <- ifelse(!is.null(doseAmt),doseAmt/AUCINF_obs,NA)
          Vz_pred <- ifelse(!is.null(doseAmt),doseAmt/(Lambda_z*AUCINF_pred),NA); Cl_pred <- ifelse(!is.null(doseAmt),doseAmt/AUCINF_pred,NA)
        }else{
          Vz_obs  <- ifelse(!is.null(doseAmt),doseAmt/(Lambda_z*AUCtau),NA);  Cl_obs  <- ifelse(!is.null(doseAmt),doseAmt/AUCtau,NA)
          Vz_pred <- NA; Cl_pred <- NA
        }
        
      }
      
      if(AUMClast != 0 & AUMCINF_obs != 0){
        AUMCINF_obs  <- AUMClast+AUMCINF_obs;  AUMCINF_D_obs  <- ifelse(!is.null(doseAmt),AUMCINF_obs/doseAmt,NA);  AUMC_pExtrap_obs  <- 100*(AUMCINF_obs-AUMClast)/AUMCINF_obs
        AUMCINF_pred <- AUMClast+AUMCINF_pred; AUMCINF_D_pred <- ifelse(!is.null(doseAmt),AUMCINF_pred/doseAmt,NA); AUMC_pExtrap_pred <- 100*(AUMCINF_pred-AUMClast)/AUMCINF_pred
      }
      
      if (AUCINF_obs  != 0 & AUMCINF_obs  != 0){MRTINF_obs  <- ifelse(adminType != "iv-infusion", AUMCINF_obs/AUCINF_obs,   ((AUMCINF_obs/AUCINF_obs)-(TI/2)))} #else{MRTINF_obs  <- "NaN"}
      if (AUCINF_pred != 0 & AUMCINF_pred != 0){MRTINF_pred <- ifelse(adminType != "iv-infusion", AUMCINF_pred/AUCINF_pred, ((AUMCINF_pred/AUCINF_pred)-(TI/2)))} #else{MRTINF_pred <- "NaN"}
      
      if(doseType == "ns" && (adminType == "iv-bolus"|adminType == "iv-infusion")){
        Vss_obs  <- MRTINF_obs*Cl_obs
        Vss_pred <- MRTINF_pred*Cl_pred
      }
      
      if (doseType == "ss" && ssnPt != 0){
        Cavg <- AUCtau/Tau
        Clss <- ifelse(!is.null(doseAmt),doseAmt/AUCtau,NA)
        Cmax <- max(nconc)
        p_Fluctuation <- 100*(Cmax-Cmin)/Cavg
        if (Lambda_z != "NaN"){Accumulation_Index <- 1/(1-exp(-Lambda_z*Tau))}
        
        # re-define MRTINF for steady state and calculate Vss
        if (AUCINF_obs  != "NaN" & AUMCINF_obs  !=  "NaN" & AUCtau != 0){
          if (adminType == "iv-infusion"){
            MRTINF_obs <- ((AUMCtau+Tau*(AUCINF_obs-AUCtau))/(AUCtau)-(TI/2)); MRTINF_pred <- ((AUMCtau+Tau*(AUCINF_pred-AUCtau))/(AUCtau)-(TI/2))
          }else{
            MRTINF_obs <- (AUMCtau+Tau*(AUCINF_obs-AUCtau))/AUCtau; MRTINF_pred <- (AUMCtau+Tau*(AUCINF_pred-AUCtau))/AUCtau
          }
          Vss_obs  <- MRTINF_obs*Clss
          Vss_pred <- MRTINF_pred*Clss
        }
      }
    }
  }
  
  if(!exists("C0")) C0 <- NA
  if(!exists("Tmax")) Tmax <- NA
  if(!exists("Cmax")) Cmax <- NA
  if(!exists("Cmax_D")) Cmax_D <- NA
  if(!exists("Tlast")) Tlast <- NA
  if(!exists("Clast")) Clast <- NA
  if(!exists("AUClast")) AUClast <- NA
  if(!exists("AUMClast")) AUMClast <- NA
  if(!exists("MRTlast")) MRTlast <- NA
  if(!exists("No_points_Lambda_z")) No_points_Lambda_z <- NA
  if(!exists("AUC_pBack_Ext_obs")) AUC_pBack_Ext_obs <- NA
  if(!exists("AUC_pBack_Ext_pred")) AUC_pBack_Ext_pred <- NA
  if(!exists("AUClower_upper")) AUClower_upper <- NA
  if(!exists("Rsq")) Rsq <- NA
  if(!exists("Rsq_adjusted")) Rsq_adjusted <- NA
  if(!exists("Corr_XY")) Corr_XY <- NA
  if(!exists("Lambda_z")) Lambda_z <- NA
  if(!exists("Lambda_z_lower")) Lambda_z_lower <- NA
  if(!exists("Lambda_z_upper")) Lambda_z_upper <- NA
  if(!exists("HL_Lambda_z")) HL_Lambda_z <- NA
  if(!exists("AUCINF_obs")) AUCINF_obs <- NA
  if(!exists("AUCINF_D_obs")) AUCINF_D_obs <- NA
  if(!exists("AUC_pExtrap_obs")) AUC_pExtrap_obs <- NA
  if(!exists("Vz_obs")) Vz_obs <- NA
  if(!exists("Cl_obs")) Cl_obs <- NA
  if(!exists("AUCINF_pred")) AUCINF_pred <- NA
  if(!exists("AUCINF_D_pred")) AUCINF_D_pred <- NA
  if(!exists("AUC_pExtrap_pred")) AUC_pExtrap_pred <- NA
  if(!exists("Vz_pred")) Vz_pred <- NA
  if(!exists("Cl_pred")) Cl_pred <- NA
  if(!exists("AUMCINF_obs")) AUMCINF_obs <- NA
  if(!exists("AUMC_pExtrap_obs")) AUMC_pExtrap_obs <- NA
  if(!exists("AUMCINF_pred")) AUMCINF_pred <- NA
  if(!exists("AUMC_pExtrap_pred")) AUMC_pExtrap_pred <- NA
  if(!exists("MRTINF_obs")) MRTINF_obs <- NA
  if(!exists("MRTINF_pred")) MRTINF_pred <- NA
  if(!exists("Tmin")) Tmin <- NA
  if(!exists("Cmin")) Cmin <- NA
  if(!exists("Cavg")) Cavg <- NA
  if(!exists("AUCtau")) AUCtau <- NA
  if(!exists("AUMCtau")) AUMCtau <- NA
  if(!exists("Clss")) Clss <- NA
  if(!exists("Vss_obs")) Vss_obs <- NA
  if(!exists("Vss_pred")) Vss_pred <- NA
  if(!exists("p_Fluctuation")) p_Fluctuation <- NA
  if(!exists("Accumulation_Index")) Accumulation_Index <- NA
  
  if (!is.null(simFile) & dset == "obs"){
    NCAprm <- as.numeric(c(C0,Tmax,0,0,0,Cmax,0,0,0,Cmax_D,Tlast,Clast,AUClast,0,0,0,AUMClast,0,0,0,MRTlast,No_points_Lambda_z,AUClower_upper,0,0,0,Rsq,Rsq_adjusted,Corr_XY,Lambda_z,Lambda_z_lower,Lambda_z_upper,HL_Lambda_z,0,0,0,AUCINF_obs,0,0,0,AUCINF_D_obs,AUC_pExtrap_obs,AUC_pBack_Ext_obs,Vz_obs,Cl_obs,AUCINF_pred,0,0,0,AUCINF_D_pred,AUC_pExtrap_pred,AUC_pBack_Ext_pred,Vz_pred,Cl_pred,AUMCINF_obs,AUMC_pExtrap_obs,AUMCINF_pred,AUMC_pExtrap_pred,MRTINF_obs,MRTINF_pred,Tau,Tmin,Cmin,Cavg,AUCtau,AUMCtau,Clss,Vss_obs,Vss_pred,p_Fluctuation,Accumulation_Index))
    
    names(NCAprm) <- c("C0","Tmax","simTmax","dTmax","npdeTmax","Cmax","simCmax","dCmax","npdeCmax","Cmax_D","Tlast","Clast","AUClast","simAUClast","dAUClast","npdeAUClast","AUMClast","simAUMClast","dAUMClast","npdeAUMClast","MRTlast","No_points_Lambda_z","AUClower_upper","simAUClower_upper","dAUClower_upper","npdeAUClower_upper","Rsq","Rsq_adjusted","Corr_XY","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","simHL_Lambda_z","dHL_Lambda_z","npdeHL_Lambda_z","AUCINF_obs","simAUCINF_obs","dAUCINF_obs","npdeAUCINF_obs","AUCINF_D_obs","AUC_pExtrap_obs","AUC_pBack_Ext_obs","Vz_obs","Cl_obs","AUCINF_pred","simAUCINF_pred","dAUCINF_pred","npdeAUCINF_pred","AUCINF_D_pred","AUC_pExtrap_pred","AUC_pBack_Ext_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred","p_Fluctuation","Accumulation_Index")
  }else{
    NCAprm <- as.numeric(c(C0,Tmax,Cmax,Cmax_D,Tlast,Clast,AUClast,AUMClast,MRTlast,No_points_Lambda_z,AUClower_upper,Rsq,Rsq_adjusted,Corr_XY,Lambda_z,Lambda_z_lower,Lambda_z_upper,HL_Lambda_z,AUCINF_obs,AUCINF_D_obs,AUC_pExtrap_obs,AUC_pBack_Ext_obs,Vz_obs,Cl_obs,AUCINF_pred,AUCINF_D_pred,AUC_pExtrap_pred,AUC_pBack_Ext_pred,Vz_pred,Cl_pred,AUMCINF_obs,AUMC_pExtrap_obs,AUMCINF_pred,AUMC_pExtrap_pred,MRTINF_obs,MRTINF_pred,Tau,Tmin,Cmin,Cavg,AUCtau,AUMCtau,Clss,Vss_obs,Vss_pred,p_Fluctuation,Accumulation_Index))
    
    names(NCAprm) <- c("C0","Tmax","Cmax","Cmax_D","Tlast","Clast","AUClast","AUMClast","MRTlast","No_points_Lambda_z","AUClower_upper","Rsq","Rsq_adjusted","Corr_XY","Lambda_z","Lambda_z_lower","Lambda_z_upper","HL_Lambda_z","AUCINF_obs","AUCINF_D_obs","AUC_pExtrap_obs","AUC_pBack_Ext_obs","Vz_obs","Cl_obs","AUCINF_pred","AUCINF_D_pred","AUC_pExtrap_pred","AUC_pBack_Ext_pred","Vz_pred","Cl_pred","AUMCINF_obs","AUMC_pExtrap_obs","AUMCINF_pred","AUMC_pExtrap_pred","MRTINF_obs","MRTINF_pred","Tau","Tmin","Cmin","Cavg","AUCtau","AUMCtau","Clss","Vss_obs","Vss_pred","p_Fluctuation","Accumulation_Index")
  }
  return(NCAprm)
}

