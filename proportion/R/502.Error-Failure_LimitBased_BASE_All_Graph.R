########################################################################################
#' Plots error, long term power and pass/fail criteria using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using 6 base methods (Wald, Wald-T, Likelihood, Score, Logit-Wald, ArcSine)
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrAll(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrAll<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errAll(n,alp,phi,f)
  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for base methods") +
    ggplot2::labs(x = "Method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity")

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using Wald method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrWD(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrWD<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errWD(n,alp,phi,f)
  errdf$method = as.factor("Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Wald method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Score method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using Score method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrSC(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrSC<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errSC(n,alp,phi,f)
  errdf$method = as.factor("Score")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Score method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using ArcSine method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using ArcSine method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrAS(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrAS<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errAS(n,alp,phi,f)
  errdf$method = as.factor("ArcSine")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for ArcSine method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Wald-T method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using Wald-T method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrTW(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrTW<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errTW(n,alp,phi,f)
  errdf$method = as.factor("Wald-T")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Wald-T method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Logit Wald method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using Logit Wald method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrLT(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrLT<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errLT(n,alp,phi,f)
  errdf$method = as.factor("Logit Wald")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Logit Wald method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}


########################################################################################
#' Plots error, long term power and pass/fail criteria using Likelihood Ratio method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @details  Plots of error, long term power and pass/fail
#' criteria using Likelihood Ratio method
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2
#' PloterrLR(n,alp,phi,f)
#' @export
##### 10. Expected length for a given n and alpha level for 6 base methods
PloterrLR<-function(n,alp,phi,f)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errLR(n,alp,phi,f)
  errdf$method = as.factor("Likelihood Ratio")

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Likelihood Ratio method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Bayesian method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @param a - Beta parameters for hypo "p"
#' @param b - Beta parameters for hypo "p"
#' @details  Plot of Bayesian Highest Probability Density
#' (HPD) and two tailed intervals using error due to the difference of achieved and
#' nominal level of significance for the \eqn{n + 1} intervals
#' for the Beta - Binomial conjugate prior model for the probability of success \code{p}
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05; phi=0.05; f=-2;a=0.5;b=0.5
#' PloterrBA(n,alp,phi,f,a,b)
#' @export
#8.BAYESIAN
PloterrBA<-function(n,alp,phi,f,a,b)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  if ((class(a) != "integer") & (class(a) != "numeric") || a<0  ) stop("'a' has to be greater than or equal to 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || b<0  ) stop("'b' has to be greater than or equal to 0")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errBA(n,alp,phi,f,a,b)

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,method=alpdf$method)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, method=thetadf$method)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = method, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Bayesian method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}

########################################################################################
#' Plots error, long term power and pass/fail criteria using Exact method
#' @param n - Number of trials
#' @param alp - Alpha value (significance level required)
#' @param phi - Null hypothesis value
#' @param f - Failure criterion
#' @param e - Exact method indicator  in [0, 1] {1: Clopper Pearson, 0.5: Mid P}
#' The input can also be a range of values between 0 and 1.
#' @details  Plots of Confidence interval for \code{p}
#' based on inverting equal-tailed binomial tests with null hypothesis \eqn{H0: p = p0}
#' using error due to the difference of achieved and nominal level of significance for the \eqn{n + 1} intervals
#' @family Error for base methods
#' @examples
#' n=20; alp=0.05;phi=0.05; f=-2;e=0.5 # Mid-p
#' PloterrEX(n,alp,phi,f,e)
#' n=20; alp=0.05;phi=0.05; f=-2;e=1 #Clopper-Pearson
#' PloterrEX(n,alp,phi,f,e)
#' n=20; alp=0.05;phi=0.05; f=-2;e=c(0.1,0.5,0.95,1) #Range including Mid-p and Clopper-Pearson
#' PloterrEX(n,alp,phi,f,e)
#' @export
#7.EXACT METHOD
PloterrEX<-function(n,alp,phi,f,e)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(alp)) stop("'alpha' is missing")
  if (missing(phi)) stop("'phi' is missing")
  if (missing(f)) stop("'f' is missing")
  if (missing(e)) stop("'e' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || n<=0 ) stop("'n' has to be greater than 0")
  if (alp>1 || alp<0 || length(alp)>1) stop("'alpha' has to be between 0 and 1")
  if (phi>1 || phi<0 || length(phi)>1) stop("Null hypothesis 'phi' has to be between 0 and 1")
  if ((class(f) != "integer") & (class(f) != "numeric")|| length(f)>1) stop("'f' has to be numeric value")
  if ((class(e) != "integer") & (class(e) != "numeric") || any(e>1) || any(e<0)) stop("'e' has to be between 0 and 1")
  if (length(e)>10) stop("'e' can have only 10 intervals")
  method=value=Fail_Pass=NULL

  #### Calling functions and creating df
  errdf=  errEX(n,alp,phi,f,e)
  errdf$e = as.factor(errdf$e)

  alpdf=  errdf[,c(1,3,4)]
  thetadf=errdf[,c(2,3,4)]
  vdfa=data.frame(value=alpdf$delalp ,mark="Increase in nominal error" ,Fail_Pass=alpdf$Fail_Pass ,e=alpdf$e)
  vdft=data.frame(value=thetadf$theta ,mark="Long term power of test",Fail_Pass=thetadf$Fail_Pass, e=thetadf$e)
  full.df=rbind(vdfa,vdft)

  ggplot2::ggplot(full.df, ggplot2::aes(x = e, y = value, fill = Fail_Pass)) +
    ggplot2::labs(title = "Error, long term power and pass/fail for Exact method") +
    ggplot2::facet_grid(mark ~ .,scales="free_y") +
    ggplot2::geom_bar(stat="identity",position = "identity",width=0.5)

}
