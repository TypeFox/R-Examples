########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted Wald-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrAWD(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrAWD<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errAWD(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted Wald method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted Score-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrASC(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED Score - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrASC<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errASC(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted Score")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted Score method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted ArcSine-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrAAS(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrAAS<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errAAS(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted ArcSine")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted ArcSine method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted Logit Wald-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrALT(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED Logit Wald - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrALT<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errALT(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted Logit Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted Logit Wald method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted Wald-T-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrATW(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrATW<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || h<0  ) stop("'h' has to be greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errATW(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted Wald-T")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted Wald-T method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using adjusted Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plot of adjusted Likelihood Ratio-type interval using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2;phi=0.99; f=-2
#' PloterrALR(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED Likelihood Ratio - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrALR<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0  || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errALR(n,alp,h,phi,f)
  errdf$method = as.factor("Adjusted Likelihood Ratio")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted Likelihood Ratio method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plot of error, long term power and pass/fail criteria using 6 adjusted methods
#' (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param h - Adding factor
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Calculates  error, long term power and pass/fail
#' criteria using 6 adjusted methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @family Error for adjusted methods
#' @examples
#' n=20; alp=0.05; h=2; phi=0.99; f=-2
#' PloterrAAll(n,alp,h,phi,f)
#' @export
##### 1.ADJUSTED Likelihood Ratio - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrAAll<-function(n,alp,h,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(h)) stop("'h' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if ((class(h) != "integer") & (class(h) != "numeric") || length(h) >1|| h<0 || !(h%%1 ==0)) stop("'h' has to be an integer greater than or equal to 0")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errAAll(n,alp,h,phi,f)
  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for adjusted methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity")

}
