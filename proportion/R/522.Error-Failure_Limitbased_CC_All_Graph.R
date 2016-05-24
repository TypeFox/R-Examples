#########################################################################################
#' Plots the calculated error, long term power and pass/fail criteria using 5 continuity
#' corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plots the  error, long term power and pass/fail criteria calculaed using 5
#' continuity corrected methods (Wald, Wald-T, Score, Logit-Wald, ArcSine)
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05; c=1/(2*n); f=-2
#' PloterrCAll(n,alp,phi,c,f)
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCAll<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCAll(n,alp,phi,c,f)
  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity")

}

################################################################################################
#' Plots the error, long term power and pass/fail criteria for continuity corrected Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plot of Wald-type interval with continuity correction using error
#' due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' PloterrCWD(n,alp,phi,c,f)
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCWD<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCWD(n,alp,phi,c,f)
  errdf$method = as.factor("Continuity corrected Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

################################################################################################
#' Plots the error, long term power and pass/fail criteria for continuity corrected ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plot of continuity corrected Wald-type interval for the arcsine
#' transformation of the parameter \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' PloterrCAS(n,alp,phi,c,f)
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCAS<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCAS(n,alp,phi,c,f)
  errdf$method = as.factor("Continuity corrected ArcSine")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected ArcSine") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

################################################################################################
#' Plots the error, long term power and pass/fail criteria for continuity corrected Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plot of continuity corrected Wald-type interval based on the logit transformation of \code{p}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' PloterrCLT(n,alp,phi,c,f)
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCLT<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCLT(n,alp,phi,c,f)
  errdf$method = as.factor("Continuity corrected Logit Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected Logit Wald") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

################################################################################################
#' Plots the error, long term power and pass/fail criteria for continuity corrected Wald-t method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plot of approximate and continuity corrected method based on a t_approximation of the standardized point estimator
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' PloterrCTW(n,alp,phi,c,f)
#' @export
##### 1.CC WALD - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCTW<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCTW(n,alp,phi,c,f)
  errdf$method = as.factor("Continuity corrected Wald-t")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected Wald-t") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

################################################################################################
#' Plots the error, long term power and pass/fail criteria for continuity corrected Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param c - Continuity correction
#' @param f - Failure criterion
#' @details  Plot of continuity corrected score test approach using error due to the
#' difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for continuity corrected methods
#' @examples
#' n=5; alp=0.05; phi=0.05;c=1/(2*n); f=-2
#' PloterrCSC(n,alp,phi,c,f)
#' @export
##### 1.CC Score - DELTA_ALPHA, THETA,F-ERROR,POWER,FAILURE
PloterrCSC<-function(n,alp,phi,c,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(c)) stop("'c' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if (c<=0 || c>(1/(2*n))) stop("'c' has to be positive and less than or equal to 1/(2*n)")
  if ((class(f) != "integer") & (class(f) != "numeric")) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errCSC(n,alp,phi,c,f)
  errdf$method = as.factor("Continuity corrected Score")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for continuity corrected Score") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}
