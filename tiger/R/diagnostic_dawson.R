diagnostic_dawson <- function(modelled, measured, p=NA, m=NA,
additional=TRUE, use_qualV=FALSE, diff.ecdf=NA){
  n = length(modelled)
  if(n!=length(measured)) stop("Modelled and measured time series must be of equal length")
  
  measured[is.na(modelled)] <- NA
  modelled[is.na(measured)] <- NA

  difff = measured - modelled
  if(typeof(diff.ecdf)!="closure"){
     diff.ecdf <- ecdf(difff)
  }
  absdifff = abs(difff)
  sqdifff = difff^2
  sumdifff = sum(difff, na.rm=TRUE)
  summeasured = sum(measured, na.rm=TRUE)
  summodelled <- sum(modelled,na.rm=TRUE)
  diff_measured <- diff(measured)
  diff_modelled <- diff(modelled)
  abs_diff_measured <- abs(measured)
  abs_diff_modelled <- abs(modelled)

  lcs_slope_f <- min(c(abs_diff_measured[abs_diff_measured>0], 
                   abs_diff_modelled[abs_diff_modelled>0]), na.rm=TRUE)
  sumabsdifff = sum(absdifff, na.rm=TRUE)
  sumsqdifff = sum(sqdifff, na.rm=TRUE)
  max_measured = max(measured,na.rm=TRUE) 
  min_measured = min(measured,na.rm=TRUE) 
  max_modelled = max(modelled, na.rm=TRUE)
  
  AME = max(absdifff, na.rm=TRUE)
  PDIFF = max_measured - max_modelled  
  MAE =  sumabsdifff/n
  ME = sumdifff / n
  RMSE = sqrt(sumsqdifff/n)
  R4MS4E = (sum(sqdifff^2, na.rm=TRUE)/n)^(1/4)
  lnRMSE=log(RMSE)
  AIC = m*lnRMSE+2*p
  BIC = m*lnRMSE+p*log(m)
  #um die anzahl der sign changes zu bestimmen:
  #vorzeichen als (0,1)
  SIGN = rep(0,n)
  SIGN[difff == 0] <- NA
  SIGN[difff > 0]<-1
  #Jeder Vorzeichenwechsel wir 1 sonst 0
  SC <- abs(diff(SIGN))
  
  NSC = sum(SC, na.rm=TRUE)
  mean_meas <- mean(measured, na.rm=TRUE)
  mean_modelled <- summodelled / n
  difff_mean <- measured - mean_meas
  difff_modelled <- modelled - mean_modelled
  absdifff_mean <- abs(difff_mean)
  sqdiff_mean <- difff_mean^2
  sqdiff_modelled <- difff_modelled^2
  sumsqdiff_mean <- sum(sqdiff_mean, na.rm=TRUE)
  sumsqdiff_modelled <- sum(sqdiff_modelled, na.rm=TRUE)
  sumabsdifff_mean <- sum(absdifff_mean, na.rm=TRUE)
  RAE <- sumabsdifff / sumabsdifff_mean
  PEP <- 100*(max_measured - max_modelled)/max_measured
  reldifff <- difff/measured
  relabsdiff <- abs(reldifff)
  MARE <- mean(relabsdiff, na.rm=TRUE)
  MdAPE <- 100*median(relabsdiff, na.rm=TRUE)
  MRE <- mean(reldifff)
  sqreldiff <- reldifff^2
  MSRE = mean(sqreldiff)
  RVE = sumdifff /summeasured
  if(sum(is.na(modelled)) |sum(is.na(measured)) | any(is.na(diff_measured)) ){
       warning("Detailed handling of NA not implemented. Returning NA
       for Rsqr, other measurse may be meaningless")
       Rsqr <- NA
       IRMSE <- NA

  } else {
       Rsqr=cor(modelled,measured)^2
       IRMSE <- RMSE/sd(diff_measured)
  }
  if(sum(!is.na(modelled)) <10 || (sd(measured, na.rm=TRUE)==0 & sd(modelled, na.rm=TRUE)==0)  ){
       t_test <- NA
  } else {
       t_test <- t.test(measured,y=modelled, paired=TRUE)$statistic
  }
  CE = 1-(sumsqdifff/sumsqdiff_mean)
  IoAd = 1- (sumsqdifff/sum((abs(modelled-mean_meas)+absdifff_mean)^2, na.rm=TRUE))
  PI = 1- (sumsqdifff/sum(diff_measured^2,na.rm=TRUE))
  MSLE = mean((log(measured) - log(modelled))^2)
  MSDE = mean((diff_measured - diff_modelled)^2)
  if(additional){
      I=summeasured/summodelled
      if(sum(is.na(modelled)) ==0  ){
         lagtim <- lagtime(modelled, y=measured)
      } else {
         lagtim <- NA
      }
      mean_reldiff <- mean(diff_measured/diff_modelled, na.rm=TRUE)
      DE <- sum( (diff_measured/diff_modelled)<0 ,na.rm=TRUE)
      krel <- k_rel(modelled, measured)
      spann <- max_measured-min_measured
      EQ <- mean(diff.ecdf(difff), na.rm=TRUE)
  }
  if(use_qualV){
         # require(qualV)
   lcs_slope <- function(o,p,f){
       time <- 1:length(o)
       o_fl <- f.slope(time,o,f=f)
       p_fl <- f.slope(time,p,f=f)
       return(LCS(o_fl,p_fl)$QSI)
   }
  }

  if(additional&&use_qualV){
      return(data.frame(AME=AME, PDIFF=PDIFF,MAE=MAE, ME=ME,RMSE=RMSE,
      R4MS4E=R4MS4E,  AIC=AIC, BIC=BIC, NSC=NSC, RAE=RAE,
      PEP=PEP,MARE=MARE, MdAPE=MdAPE,MRE=MRE,MSRE=MSRE,RVE=RVE, Rsqr=Rsqr,
      CE=CE, IoAd=IoAd, PI=PI, MSLE=MSLE, MSDE=MSDE, t_test=t_test,
      IRMSE=IRMSE, I=I, lagtime=lagtim, rel_diff=mean_reldiff, DE=DE,
      krel=krel, span=spann,
      CMAE=CMAE(o=modelled, p=measured),
    CMSE=CMSE(o=measured, p=modelled),
    MAGE=MAGE(o=measured, p=modelled),
    MALE=MALE(o=measured, p=modelled),
    MAOE=MAOE(o=measured, p=modelled),
    MAPE=MAPE(o=measured, p=modelled),
    MSE=MSE(o=measured, p=modelled),
    MSOE=MSOE(o=measured, p=modelled),
    RCMSE=RCMSE(o=measured, p=modelled),
    lcs_slope=lcs_slope(o=measured, p=modelled, lcs_slope_f),
    RMSGE=RMSGE(o=measured, p=modelled),
    RMSLE=RMSLE(o=measured, p=modelled),
    RMSOE=RMSOE(o=measured, p=modelled),
    RSMSE=if(sum(!is.na(measured))==0){0}else{RSMSE(o=measured, p=modelled)},
    RSMSGE=RSMSGE(o=measured, p=modelled),
    RSMSLE=RSMSLE(o=measured, p=modelled),
    SMAE=if(sum(!is.na(measured))==0){0}else{SMAE(o=measured, p=modelled)},
    SMAGE=SMAGE(o=measured, p=modelled),
    SMALE=SMALE(o=measured, p=modelled),
    SMSE=if(sum(!is.na(measured))==0){0}else{SMSE(o=measured, p=modelled)},
    SMSLE=SMSLE(o=measured, p=modelled),
    GRI=GRI(o=measured, p=modelled), EQ=EQ
      ))
  } else if(additional) {
      return(data.frame(AME=AME, PDIFF=PDIFF,MAE=MAE, ME=ME,RMSE=RMSE,
      R4MS4E=R4MS4E,  AIC=AIC, BIC=BIC, NSC=NSC, RAE=RAE,
      PEP=PEP,MARE=MARE, MdAPE=MdAPE,MRE=MRE,MSRE=MSRE,RVE=RVE, Rsqr=Rsqr,
      CE=CE, IoAd=IoAd, PI=PI, MSLE=MSLE, MSDE=MSDE, t_test=t_test,
      IRMSE=IRMSE, I=I, lagtime=lagtim, rel_diff=mean_reldiff, DE=DE,
      krel=krel, span=spann, EQ=EQ))
  } else if(use_qualV) {
      return(data.frame(AME=AME, PDIFF=PDIFF,MAE=MAE, ME=ME,RMSE=RMSE,
      R4MS4E=R4MS4E,  AIC=AIC, BIC=BIC, NSC=NSC, RAE=RAE,
      PEP=PEP,MARE=MARE, MdAPE=MdAPE,MRE=MRE,MSRE=MSRE,RVE=RVE, Rsqr=Rsqr,
      CE=CE, IoAd=IoAd, PI=PI, MSLE=MSLE, MSDE=MSDE, t_test=t_test,
      IRMSE=IRMSE, 
      CMAE=CMAE(o=modelled, p=measured),
    CMSE=CMSE(o=measured, p=modelled),
    MAGE=MAGE(o=measured, p=modelled),
    MALE=MALE(o=measured, p=modelled),
    MAOE=MAOE(o=measured, p=modelled),
    MAPE=MAPE(o=measured, p=modelled),
    MSE=MSE(o=measured, p=modelled),
    MSOE=MSOE(o=measured, p=modelled),
    RCMSE=RCMSE(o=measured, p=modelled),
    lcs_slope=lcs_slope(o=measured, p=modelled, lcs_slope_f),
    RMSGE=RMSGE(o=measured, p=modelled),
    RMSLE=RMSLE(o=measured, p=modelled),
    RMSOE=RMSOE(o=measured, p=modelled),
    RSMSE=if(sum(!is.na(measured))==0){0}else{RSMSE(o=measured, p=modelled)},
    RSMSGE=RSMSGE(o=measured, p=modelled),
    RSMSLE=RSMSLE(o=measured, p=modelled),
    SMAE=if(sum(!is.na(measured))==0){0}else{SMAE(o=measured, p=modelled)},
    SMAGE=SMAGE(o=measured, p=modelled),
    SMALE=SMALE(o=measured, p=modelled),
    SMSE=if(sum(!is.na(measured))==0){0}else{SMSE(o=measured, p=modelled)},
    SMSLE=SMSLE(o=measured, p=modelled),
    GRI=GRI(o=measured, p=modelled)
      ))
  } else {
      return(data.frame(AME=AME, PDIFF=PDIFF,MAE=MAE, ME=ME,RMSE=RMSE,
      R4MS4E=R4MS4E,  AIC=AIC, BIC=BIC, NSC=NSC, RAE=RAE,
      PEP=PEP,MARE=MARE, MdAPE=MdAPE,MRE=MRE,MSRE=MSRE,RVE=RVE, Rsqr=Rsqr,
      CE=CE, IoAd=IoAd, PI=PI, MSLE=MSLE, MSDE=MSDE, t_test=t_test,
      IRMSE=IRMSE))
  }
  
}
