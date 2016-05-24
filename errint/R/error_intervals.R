#' Error Intervals
#'
#' @rdname error_interval
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{error_interval} creates an object of class
#'  \code{c("error_interval", "list")}.
#'
#' @param phi a vector with residual values used to compute the error interval.
#' @param s confidence level, e,g. s=0.05 for the standard 95 percent confidence interval.
#' @param dist assumed distribution for the noise in the data.
#' @param tol used to normalize residual values to (0,1) when beta is the
#' assumed distribution. The formula used is abs(phi)/(max(abs(phi))+tol).
#' @param ... additional arguments to be passed to the low level
#' error_interval building functions (see below).
#'
#' @return \code{error_interval} returns an object of class \code{c("error_interval","list")}
#' with information regarding the error intervals built.
#' @export
#'
#' @examples
#' error_interval(rnorm(100))
#'
#' error_interval(rnorm(100),s=0.1,dist="lm")


error_interval <- function(phi, s=0.05, dist="n",tol=10^-6, ...) UseMethod("error_interval");

#' Error Intervals
#'
#' @rdname error_interval
#'
#' @description \code{as.error_interval} attempts to coerce its argument \code{x} into an
#' object of class \code{c("error_interval", class(x))}. If this is not possible
#' \code{x} is returned unchanged.
#'
#' @param x an R object.
#'
#' @return \code{as.error_interval} returns an object of class
#' \code{c("error_interval",class(x))} with information contained in x if possible.
#' Returns x otherwise.
#' @export
#'
#' @examples
#'
#'
#' l<-list(min=-1,max=1,err=0.05,s=0.1,dist="n",phi=rnorm(1000))
#' as.error_interval(l)
#'
#' v<-c("a","b")
#' as.error_interval(v)
as.error_interval <- function (x){
  if (is.error_interval(x)){
    return(x);
  }
  if (is.list(x)){
    res<-x;
    class(res)<-c("error_interval",class(x));
    return(res);
  }
  else{
    warning("It is not possible to coerce x to class error_interval. x returned unchanged");
    return(x);
  }
}

#' Error Intervals
#'
#' @rdname error_interval
#'
#' @description \code{is.error_interval} returns TRUE if x is an R object with
#' \code{"error_interval"} as one of its classes. It returns FALSE otherwise.
#'
#' @return  \code{is.error_interval}returns TRUE if x is an R object with
#' \code{"error_interval"} as one of its classes. FALSE otherwise.
#' @export
#'
#' @examples
#'
#'
#' l<-list(min=-1,max=1,err=0.05,s=0.1,dist="n",phi=rnorm(1000))
#' is.error_interval(l)
#' res<-as.error_interval(l)
#' is.error_interval(res)
is.error_interval <- function (x){
  return("error_interval" %in% class(x));
}


#' Error Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{error_interval.default} creates an object of class
#'  \code{c("error_interval", "list")}.
#'
#' @param phi a vector with residual values used to compute the error interval.
#' @param s confidence level, e,g. s=0.05 for the standard 95 percent confidence interval.
#' @param dist assumed distribution for the noise in the data.
#' @param tol used to normalize residual values to (0,1) when beta is the
#' assumed distribution. The formula used is abs(phi)/(max(abs(phi))+tol).
#' @param ... additional arguments to be passed to the low level
#' error_interval building functions (see below).
#'
#' @return Returns an object of class \code{c("error_interval","list")}
#' with information regarding the error intervals built.
#' @export
#'
#' @examples
#' error_interval(rnorm(100))
#'
#' error_interval(rnorm(100),s=0.1,dist="lm")


error_interval.default <- function(phi, s=0.05, dist="n",tol=10^-6, ...)
{
  res<- if (dist=="l") {int_lap(phi,s,...);}
  else if (dist=="lm") {int_lap_mu(phi,s,...);}
  else if (dist=="n") {int_gau(phi,s,...);}
  else if (dist=="nm") {int_gau_mu(phi,s,...);}
  else if (dist=="b") {int_beta(abs(phi)/(max(abs(phi))+tol),s,...);}
  else if (dist=="w") {int_weibull(abs(phi),s,...);}
  return(res);
}

#' Printing Error Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{print} objects of class \code{error_interval}.
#'
#' @param x object of class \code{error_interval} to be printed.
#' @param ... optional arguments.
#'
#' @seealso \link{error_interval}
#' @export
#'
#' @examples
#' res<-error_interval(rnorm(100))
#' print(res)


print.error_interval<- function(x, ...){
  cat(sprintf("(%s,%s)\n",round(x[["min"]], digits=4),round(x[["max"]], digits=4)));
}

#' Error Intervals Summaries
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{summary} produces summaries for
#' objects of class \code{error_interval}.
#'
#' @param object object of class \code{error_interval} to be printed.
#' @param ... optional arguments.
#'
#' @return Object of class \code{c("summary.error_interval","list")} corresponding
#' to the summary of x.
#' @seealso \link{error_interval}
#' 
#' @S3method summary error_interval
#'
#' @examples
#' res<-error_interval(rnorm(100))
#' summary(res)

summary.error_interval<- function(object, ...){
  range<-object[["max"]]-object[["min"]];
  rel_err<-object[["err"]]/(1-object[["s"]]);
  dist<-if (object[["dist"]]=="l") {dist<-"Zero-mu Laplace";}
  else if (object[["dist"]]=="lm") {dist<-"General Laplace";}
  else if (object[["dist"]]=="n") {dist<-"Zero-mu Gaussian";}
  else if (object[["dist"]]=="nm") {dist<-"General Gaussian";}
  else if (object[["dist"]]=="b") {dist<-"Beta";}
  else if (object[["dist"]]=="w") {dist<-"Weibull";}
  res<-list(min=object[["min"]],max=object[["max"]],err=object[["err"]],
            s=object[["s"]],dist=dist,range=range,rel_err=rel_err);
  class(res)<-c("summary.error_interval","list");
  return(res);
}

#' Printing Error Intervals Summaries
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{print} objects of class \code{summary.error_interval}.
#'
#' @param x object of class \code{summary.error_interval} to be printed.
#' @param ... optional arguments.
#'
#' @seealso \link{summary} \link{error_interval}
#' @export
#'
#' @examples
#' res<-error_interval(rnorm(100))
#' summary(res)

print.summary.error_interval<- function(x, ...){
  cat("Assumed distribution:\n");
  cat(x[["dist"]]);
  cat("\n\n");
  cat("Error interval:\n");
  cat(sprintf("(%s,%s)\n\n",round(x[["min"]], digits=4),round(x[["max"]], digits=4)));
  cat("Interval range:\n");
  cat(x[["range"]]);
  cat("\n\n");
  cat("Absolute error in integral equation:\n");
  cat(sprintf("%e",x[["err"]]));
  cat("\n\n");
  cat(sprintf("Relative error in integral equation (1-s=%s):\n",1-x[["s"]]));
  cat(sprintf("%e\n",x[["rel_err"]]));
}



#' Measures
#'
#' @rdname measure
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{measure} creates an object of class
#'  \code{c("measure", "list")}.
#'
#' @param s confidence level, e,g. s=0.05 for the standard 95 percent confidence interval.
#' @param acc accuracy achieved by error intervals.
#' @param f function used to compute error of intervals. See also 'Details'.
#'
#' @return \code{measure} returns an object of class \code{c("measure","list")}
#' with information regarding the error of a set of intervals.
#' @export
#'
#' @examples
#' measure(0.1,0.7)
#'
#' measure(0.1,0.7,function(x,y){y-x})
measure <- function(s,acc,f=function(x,y){abs(x-y)}) UseMethod("measure");

#' Measures
#'
#' @rdname measure
#'
#' @description \code{as.measure} attempts to coerce its argument \code{x} into an
#' object of class \code{c("measure", class(x))}. If this is not possible
#' \code{x} is returned unchanged.
#'
#' @param x an R object.
#'
#' @return \code{as.measure} returns an object of class
#' \code{c("measure",class(x))} with information contained in x if possible.
#' Returns x otherwise.
#' @export
#'
#' @examples
#'
#'
#' l<-list(s=0.1,acc=0.78,f=function(x,y){abs(x-y)},err=0.02)
#' as.measure(l)
#'
#' v<-c("a","b")
#' as.measure(v)
as.measure <- function (x){
  if (is.measure(x)){
    return(x);
  }
  else if (is.list(x)){
    res<-x;
    class(res)<-c("measure",class(x));
    return(res);
  }
  else{
    warning("It is not possible to coerce x to class measure. x returned unchanged");
    return(x);
  }
}

#' Measures
#'
#' @rdname measure
#'
#' @description \code{is.measure} returns TRUE if x is an R object with
#' \code{"measure"} as one of its classes. It returns FALSE otherwise.
#'
#' @return  \code{is.measure} returns TRUE if x is an R object with
#' \code{"measure"} as one of its classes. FALSE otherwise.
#' @export
#'
#' @examples
#'
#'
#' l<-list(s=0.1,acc=0.78,f=function(x,y){abs(x-y)},err=0.02)
#' is.measure(l)
#' res<-as.measure(l)
#' is.measure(res)
is.measure <- function (x){
  return("measure" %in% class(x));
}

#' Measure
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{measure} creates an object of class
#'  \code{c("measure", "list")}.
#'
#' @param s confidence level, e,g. s=0.05 for the standard 95 percent confidence interval.
#' @param acc accuracy achieved by error intervals.
#' @param f function used to compute error of intervals. See also 'Details'.
#'
#' @return Returns an object of class \code{c("measure","list")}
#' with information regarding the error of a set of intervals.
#' @export
#'
#' @examples
#' measure(0.1,0.7)
#'
#' measure(0.1,0.7,function(x,y){y-x})
measure.default <- function(s,acc,f=function(x,y){abs(x-y)})
{
  res<-list(s=s,acc=acc,f=f,err=f(1-2*s,acc))
  class(res)<-c("measure","list");
  return(res);
}

#' Printing Measures
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{print} objects of class \code{measure}.
#'
#' @param x object of class \code{measure} to be printed.
#' @param ... optional arguments.
#'
#' @seealso \link{measure}
#' @export
#'
#' @examples
#' res<-measure(0.1,0.7)
#' print(res)
print.measure<- function(x, ...){
  cat(x[["err"]]*100);
  cat("%\n");
}

#' Measures Summaries
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{summary} produces summaries for
#' objects of class \code{measure}.
#'
#' @param object object of class \code{measure} to be printed.
#' @param ... optional arguments.
#'
#' @return Object of class \code{c("summary.measure","list")} corresponding
#' to the summary of x.
#' @seealso \link{measure}
#' 
#' @S3method summary measure
#
#' @examples
#' res<-measure(0.1,0.7)
#' summary(res)
summary.measure<- function(object, ...){
  res<-list(goal=1-2*object[["s"]],achieved=object[["acc"]],err=object[["err"]],
            f=object[["f"]]);
  class(res)<-c("summary.measure","list");
  return(res);
}

#' Printing Measures Summaries
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{print} objects of class \code{summary.measure}.
#'
#' @param x object of class \code{summary.measure} to be printed.
#' @param ... optional arguments.
#'
#' @seealso \link{summary} \link{measure}
#' @export
#'
#' @examples
#' res<-measure(0.1,0.7)
#' summary(res)
print.summary.measure<- function(x, ...){
  cat("Goal accuracy:\n");
  cat(x[["goal"]]*100);
  cat("%\n\n");
  cat("Achieved accuracy:\n");
  cat(x[["achieved"]]*100);
  cat("%\n\n");
  cat("Error function:\n");
  print(x[["f"]])
  cat("\n");
  cat("Error:\n");
  cat(x[["err"]]*100);
  cat("%\n");
}

#' Data Frames of Intervals
#'
#' @rdname df_intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{df_intervals} creates an object of class
#'  \code{c("df_intervals", "data.frame")}.
#'
#' @param distributions vector containing the names of the distribution
#' correspondind to each error.
#' @param errs vector of errors associated to intervals built under
#' a particular distribution assumption indicated by 'distributions'.
#'
#' @return \code{df_intervals} returns an object of class
#' \code{c("df_intervals", "data.frame")} with information regarding the error of
#' intervals built under different distribution assumptions.
#' @export
#'
#' @examples
#' df_intervals("l",0.1)
#'
#' df_intervals(c("l","lm","n","nm","b","w"),rep(0.1,6))
df_intervals <- function(distributions,errs) UseMethod("df_intervals");

#' Data Frames of Intervals
#'
#' @rdname df_intervals
#'
#' @description \code{as.df_intervals} attempts to coerce its argument \code{x} into an
#' object of class \code{c("df_intervals", class(x))}. If this is not possible
#' \code{x} is returned unchanged.
#'
#' @param x an R object.
#'
#' @return \code{as.df_intervals} returns an object of class
#' \code{c("df_intervals",class(x))} with information contained in x if possible.
#' Returns x otherwise.
#' @export
#'
#' @examples
#'
#'
#' df<-data.frame(distribution=rnorm(1000),error=rnorm(1000))
#' as.df_intervals(df)
#'
#' v<-c("a","b")
#' as.df_intervals(v)
as.df_intervals <- function (x){
  if (is.df_intervals(x)){
    return(x);
  }
  if (is.data.frame(x)){
    res<-x;
    class(res)<-c("df_intervals",class(res));
    return(res);
  }
  else{
    warning("It is not possible to coerce x to class df_intervals. x returned unchanged");
    return(x);
  }
}

#' Data Frames of Intervals
#'
#' @rdname df_intervals
#'
#' @description \code{is.df_intervals} returns TRUE if x is an R object with
#' \code{"df_intervals"} as one of its classes. It returns FALSE otherwise.
#'
#' @return  \code{is.df_intervals} returns TRUE if x is an R object with
#' \code{"df_intervals"} as one of its classes. FALSE otherwise.
#' @export
#'
#' @examples
#'
#'
#' df<-data.frame(distribution=rnorm(1000),error=rnorm(1000))
#' is.df_intervals(df)
#' res<-as.df_intervals(df)
#' is.df_intervals(res)
is.df_intervals<- function (x){
  return("df_intervals" %in% class(x));
}

#' Data Frames of Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{df_intervals} creates an object of class
#'  \code{c("df_intervals", "data.frame")}.
#'
#' @param distributions vector containing the names of the distribution
#' correspondind to each error.
#' @param errs vector of errors associated to intervals built under
#' a particular distribution assumption indicated by 'distributions'.
#'
#' @return Returns an object of class \code{c("df_intervals", "data.frame")}
#' with information regarding the error of intervals built under
#' different distribution assumptions.
#' @export
#'
#' @examples
#' df_intervals("l",0.1)
#'
#' df_intervals(c("l","lm","n","nm","b","w"),rep(0.1,6))
df_intervals.default <- function(distributions,errs){
  res<-data.frame(distribution=distributions,error=errs);
  class(res)<-c("df_intervals","data.frame");
  return(res);
}



#' Printing Data Frames of Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{print} objects of class \code{df_interval}.
#'
#' @param x object of class \code{df_interval} to be printed.
#' @param ... optional arguments.
#'
#' @seealso \link{df_intervals}
#' @export
#'
#' @examples
#' res<-df_intervals(c("l","lm","n","nm","b","w"),rep(0.1,6))
#' print(res)
print.df_intervals<- function(x, ...){
  x$error<-paste(round(x$error*100,digits=2),"%",sep="");
  print.data.frame(x);
}


#' Probability Density Functions
#'
#' @rdname prob_dens_func
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{p_laplace} computes the probability density function
#' of a random variable that has a Laplace distribution with parameters \eqn{\mu} and
#' \eqn{\sigma}.
#'
#' @param x vector of points which values we want to compute.
#' @param mu location or mean parameter of the Laplace or Gaussian distribution,
#' respectively.
#' @param sigma scale parameter of the Laplace distribution.
#'
#' @return Returns a \code{numeric} object corresponding to the value
#' of the probability density function for the given x and distribution parameters.
#'
#' @seealso \link{dlaplace}
#' @export
#' 
#' @import VGAM
#'
#' @examples
#' p_laplace(0.3)
#' p_laplace(0.3,mu=0.35,sigma=0.2)

p_laplace <- function(x,mu=0,sigma=1){
  return(VGAM::dlaplace(x,mu,sigma));
}

#' Gaussian Probability Density Function
#'
#' @rdname prob_dens_func
#'
#' @description \code{p_gaussian} computes the probability density function
#' of a random variable that has a Gaussian distribution with parameters \eqn{\mu} and
#' \eqn{{\sigma}^2}{\sigma^2}.
#'
#' @param sigma_cuad variance parameter of the Gaussian distribution.
#'
#' @seealso \link{dnorm}
#' @export
#'
#' @examples
#'
#'
#' p_gaussian(0.3)
#' p_gaussian(0.3,mu=0.35,sigma_cuad=0.2)
p_gaussian<- function(x,mu=0,sigma_cuad=1){
  return(dnorm(x,mu,sqrt(sigma_cuad)));
}

#' Beta Probability Density Function
#'
#' @rdname prob_dens_func
#'
#' @description \code{p_beta} computes the probability density function
#' of a random variable that has a Beta distribution with parameters \eqn{\alpha} and
#' \eqn{\beta}.
#'
#' @param alpha shape1 parameter of the Beta distribution.
#' @param beta shape2 parameter of the Beta distribution.
#'
#' @seealso \link{dbeta}
#' @export
#'
#' @examples
#'
#'
#' p_beta(0.3)
#' p_beta(0.3,alpha=0.35,beta=0.2)
p_beta<- function(x,alpha=1,beta=1){
  return(dbeta(x,alpha,beta));
}

#' Weibull Probability Density Function
#'
#' @rdname prob_dens_func
#'
#' @description \code{p_weibull} computes the probability density function
#' of a random variable that has a Weibull distribution with parameters \eqn{\kappa}
#' and \eqn{\lambda}.
#'
#' @param k shape parameter of the Weibull distribution.
#' @param lambda scale parameter of the Weibull distribution.
#'
#' @seealso \link{dweibull}
#' @export
#'
#' @examples
#'
#'
#' p_weibull(0.3)
#' p_weibull(0.3,k=0.35,lambda=0.2)
p_weibull<- function(x,k=1,lambda=1){
  return(dweibull(x,k,lambda));
}



#' Building Error Intervals
#'
#' @rdname build_err_int
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{int_lap} computes the error interval of a set of residuals
#' assuming a Laplace distribution with zero location for the noise.
#'
#' @param phi residual values used to compute the error interval.
#' @param s confidence level, e,g. s=0.05 for the standard 95 percent confidence interval.
#'
#' @return Returns an object of class \code{c("error_interval","list")}
#' with information of the corresponding error interval.
#'
#' @details For the Zero-\eqn{\mu} Laplace distribution the value of the corresponding integral
#' equation has a closed solution of the form \eqn{ps=-\sigma \log{2s}}{ps=-\sigma log2s}.
#'
#' @seealso \link{error_interval}
#'
#' \link{p_laplace}
#' @export
#'
#' @examples
#' int_lap(rnorm(100),0.1)
#' int_lap(rbeta(100,0.1,0.2),0.6)
int_lap <- function(phi,s){
  sigma<-mean(abs(phi),na.rm=T);
  ps<-(-sigma*log(2*s));
  res<-list(min=-ps,max=ps,err=0,s=s,dist="l",phi=phi);
  class(res)<-c("error_interval","list");
  return(res);
}

#' Zero-mu Gaussian Error Interval
#'
#' @rdname build_err_int
#'
#' @description \code{int_gau} computes the error interval of a set of residuals
#' assuming a Gaussian distribution with zero mean for the noise.
#'
#' @param ps minimum value to search for solution of the integral equation to solve.
#' See also 'Details'.
#' @param threshold step size to increase ps after each iterarion. See also 'Details'.
#' @param upper maximum value to search for solution of the integral equation to solve.
#' See also 'Details'.
#'
#' @details For the other distributions, starting with the initial value of \code{ps}
#' passed as argument, the value, \code{integral}, of the corresponding integral expression is
#' computed (see also 'References' for an in-depth explanation of this integral expression).
#'  If \code{integral} is smaller than \code{1-s} then \code{ps} is increased
#' by a step size of \code{threshold} value and \code{integral} is recomputed.
#' If \code{integral} is greater or equal than 0 or if \code{ps} gets bigger than
#' \code{upper}, the loop stops and the last value of \code{ps} will be its final value.
#'
#' @seealso \link{p_gaussian}
#' @export
#'
#' @examples
#'
#'
#' int_gau(rnorm(100),0.1)
#' int_gau(rnorm(100),0.1,0.1,10^-3,10^2)
int_gau <- function(phi,s,ps=0,threshold=10^-2,upper=10^6){
  sigma<-mean(phi^2,na.rm=T);

  while(ps < upper){
    integral<-integrate(p_gaussian,-Inf,ps,sigma_cuad=sigma,mu=0)$value;
    diff<-integral - (1-s);

    if (diff > 0){
      break;
    }
    else{
      ps<-ps+threshold;
    }
  }
  res<-list(min=-ps,max=ps,err=diff,s=s,dist="n",phi=phi);
  class(res)<-c("error_interval","list");
  return(res);
}

#' General Laplace Error Interval
#'
#' @rdname build_err_int
#'
#' @description \code{int_lap_mu} computes the error interval of a set of residuals
#' assuming a Laplace distribution.
#'
#' @export
#'
#' @examples
#'
#'
#' int_lap_mu(rnorm(100),0.1)
#' int_lap_mu(rnorm(100),0.1,0.1,10^-3,10^2)
int_lap_mu <- function(phi,s,ps=median(phi,na.rm=T),threshold=10^-2,upper=10^6){
  mu<-median(phi,na.rm=T)
  sigma<-mean(abs(phi-mu),na.rm=T);

  while(ps < upper){
    integral<-integrate(p_laplace,-Inf,ps,sigma=sigma,mu=mu)$value;
    diff<-integral - (1-s);

    if (diff > 0){
      break;
    }
    else{
      ps<-ps+threshold;
    }
  }
  res<-list(min=(2*mu-ps),max=ps,err=diff,s=s,dist="lm",phi=phi);
  class(res)<-c("error_interval","list");
  return(res);
}

#' General Gaussian Error Interval
#'
#' @rdname build_err_int
#'
#' @description \code{int_gau_mu} computes the error interval of a set of residuals
#' assuming a Gaussian distribution.
#'
#' @export
#'
#' @examples
#'
#'
#' int_gau_mu(rnorm(100),0.1)
#' int_gau_mu(rnorm(100),0.1,0.1,10^-3,10^2)
int_gau_mu <- function(phi,s,ps=mean(phi,na.rm=T),threshold=10^-2,upper=10^6){
  mu<-mean(phi,na.rm=T)
  sigma<-mean((phi-mu)^2,na.rm=T);

  while(ps < upper){
    integral<-integrate(p_gaussian,-Inf,ps,sigma_cuad=sigma,mu=mu)$value;
    diff<-integral - (1-s);

    if (diff >= 0){
      break;
    }
    else{
      ps<-ps+threshold;
    }
  }
  res<-list(min=(2*mu-ps),max=ps,err=diff,s=s,dist="nm",phi=phi);
  class(res)<-c("error_interval","list");
  return(res);
}

#' Beta Error Interval
#'
#' @rdname build_err_int
#'
#' @description \code{int_beta} computes the error interval of a set of residuals
#' assuming a Beta distribution.
#'
#' @param m1 first moment of the residuals. Used to compute \code{alpha_0}.
#' @param m2 second moment of the residuals. Used to compute \code{beta_0}.
#' @param alpha_0 initial value for Newton-Raphson method for the parameter \eqn{\alpha}.
#' See also 'Details' and \link{multiroot}.
#' @param beta_0 initial value for Newton-Raphson method for the parameter \eqn{\beta}.
#' See also 'Details' and \link{multiroot}.
#'
#' @details In addition, for the Beta distribution values of parameters \eqn{\alpha} and
#' \eqn{\beta} are estimated using Newton-Raphson method, and for the Weibull distribution
#' value of parameter \eqn{\kappa} is estimated using Newton-Raphson method and then estimated
#' value of \eqn{\lambda} is computed using a closed form that depends on \eqn{\kappa}.
#'
#' See also 'References'.
#'
#' @seealso \link{p_beta}
#' @export
#' 
#' @import rootSolve
#'
#' @examples
#'
#'
#' int_beta(runif(100,0,0.99),0.1)
#' int_beta(runif(100,0,0.99),0.1,alpha_0=1,beta_0=1)
int_beta <- function(phi,s,ps=10^-4,threshold=10^-4,upper=1,
                     m1=mean(phi,na.rm=T),
                     m2=mean(phi^2,na.rm=T),
                     alpha_0=(m1*(m1-m2))/(m2 - m1^2),
                     beta_0=(alpha_0*(1-m1)/m1)){


  c1<-mean(log(phi),na.rm=T);
  c2<-mean(log(1-phi),na.rm=T);

  f <- function(x){
    return(c(digamma(x[1])-digamma(x[1]+x[2])-c1,digamma(x[2])-digamma(x[1]+x[2])-c2));
  }


  roots<-rootSolve::multiroot(f,c(alpha_0,beta_0))$root;
  alpha<-roots[1];
  beta<-roots[2];

  while(ps < upper){
    integral<-integrate(p_beta,0,ps,alpha=alpha,beta=beta)$value;
    diff<-integral - (1-2*s);

    if (diff > 0){
      break;
    }
    else{
      ps<-ps+threshold;
    }
  }
  res<-list(min=0,max=ps,err=diff,s=s,dist="b",phi=phi);
  class(res)<-c("error_interval","list");
  return(res);
}

#' Weibull Error Interval
#'
#' @rdname build_err_int
#'
#' @description \code{int_weibull} computes the error interval of a set of residuals
#' assuming a Weibull distribution.
#'
#' See also 'Details'.
#'
#' @param k_0 initial value for Newton-Raphson method for the parameter \eqn{\kappa}.
#' See also 'Details' and \link{multiroot}.
#'
#' @seealso \link{p_weibull}
#'
#' \link{multiroot}
#' @export
#'
#' @examples
#'
#'
#' int_weibull(abs(rnorm(100)),0.1)
#' int_weibull(abs(rnorm(100)),0.1,k_0=2)
int_weibull<- function(phi,s,ps=10^-4,threshold=10^-2,upper=10^6,k_0=1){
  c1<-mean(log(phi),na.rm=T);

  f <- function(x){
    return(sum((phi^x)*log(phi))/sum(phi^x)-(1/x)-c1);
  }



  roots<-rootSolve::multiroot(f,k_0)$root;
  k<-roots[1];
  lambda<-(mean(phi^k,na.rm=T))^(1/k);

  while(ps < upper){
    integral<-integrate(p_weibull,0,ps,k=k,lambda=lambda)$value;
    diff<-integral - (1-2*s);

    if (diff > 0){
      break;
    }
    else{
      ps<-ps+threshold;
    }
  }
  res<-list(min=0,max=ps,err=diff,s=s,dist="w",phi=phi)
  class(res)<-c("error_interval","list");
  return(res);
}

#' Accuracy of Error intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{int_intervals} computes the real accuracy of a given error intervals
#' for a particular set of errors and a particular error function.
#'
#' @param interv error interval.
#' @param errors set of errors.
#' @param f error function to be used to compute error between real \code{x}
#' (\code{interv}) and predicted \code{y} (\code{errors}) values. See also 'Details'.
#' @param tol used to normalize residual values to (0,1) when beta is the
#' assumed distribution. See also 'Details'.
#'
#' @return Returns an object of class \code{c("measure","list")} with information of
#' the interval accuracy.
#'
#' @details \code{f} must be a function that takes two arguments, \code{x} and \code{y},
#' and return a numeric value.
#'
#' The formula used to normalize residual values to (0,1) when a Beta distribution is
#' assumed is \eqn{\frac{|\phi|}{max{|\phi|}+tol}}{|\phi|*(max(|\phi|)+tol)^(-1)}.
#'
#' @seealso \link{measure} \link{error_interval}
#' @export
#'
#' @examples
#' interv<-int_gau(rnorm(1000),0.1)
#' acc_intervals(interv,rnorm(1000))
#' acc_intervals(interv,rnorm(1000),function(x,y){x-y})
acc_intervals <- function(interv,errors,f=function(x,y){abs(x-y)},tol=10^-8){
  s<-interv[["s"]];
  dist<-interv[["dist"]];
  phi<-interv[["phi"]];

  acc<-if ((dist=="l") || (dist=="lm") || (dist=="n")  || (dist=="nm")){
    sum(errors>= interv$min & errors <= interv$max)/length(errors);
  }

  else if (dist=="b"){
    err_temp<-abs(errors)/(max(abs(phi))+tol);
    sum(err_temp>= interv$min & err_temp <= interv$max)/length(err_temp);
  }

  else if (dist=="w"){
    err_temp<-abs(errors);
    sum(err_temp>= interv$min & err_temp <= interv$max)/length(err_temp);
  }

  else {
    print(sprintf("Distribution %s is not a correct type of distribution",dist));
    return();
  }


  return(measure(s,acc,f));
}

#' Distribution with Best Error Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{best_distribution} computes the distribution assumption that
#' gives error intervals with the lower accuracy error for a given set of residuals.
#'
#' @param phi residual values used to compute the error interval.
#' @param errors set of real errors corresponding to the predictions of a particular
#' model.
#' @param dists character vector with the distribution assumptions to test. See also
#' 'Details'.
#' @param ... additional arguments to be passed to functions \code{error_interval}
#' and \code{acc_intervals}.
#'
#' @return Returns an object of class \code{c("df_intervals", "data.frame")} with
#' information of the distribution assumption with lower accuracy error.
#'
#' @details Allowed distribution assumptions are:
#' \itemize{
#'  \item{"n": }{Zero-mu Gaussian}
#'  \item{"nm": }{General Gaussian}
#'  \item{"l": }{Zero-mu Laplace}
#'  \item{"lm": }{General Laplace}
#'  \item{"b": }{Beta}
#'  \item{"w": }{Weibull}
#' }
#'
#' @seealso \link{df_intervals} \link{error_interval} \link{acc_intervals}
#' @export
#'
#' @examples
#' best_distribution(rnorm(10000),rnorm(10000),dists=c("n","b"))
#' best_distribution(rnorm(10000),rnorm(10000))

best_distribution <- function(phi,errors,dists=c("n","nm","l","lm","w","b"),...){
  distributions<-dists;
  errs<-vector();
  for (distribution in distributions){
    interv<-error_interval(phi, dist=distribution, ...);
    errs<-c(errs,acc_intervals(interv,errors,...)$err);
  }
  best_index<-which.min(errs);
  res<-df_intervals(distributions[best_index],errs[best_index]);
  return(res);
}

#' Sort Distributions by Better Error Intervals
#'
#' @author Jesus Prada, \email{jesus.prada@@estudiante.uam.es}
#'
#' @references 
#' Link to the scientific paper
#' 
#'Prada, Jesus, and Jose Ramon Dorronsoro. "SVRs and Uncertainty Estimates in Wind 
#'Energy Prediction." Advances in Computational Intelligence. Springer International 
#'Publishing, 2015. 564-577,
#'
#'with theoretical background for this package is provided below.
#' 
#' \url{http://link.springer.com/chapter/10.1007/978-3-319-19222-2_47}
#'
#' @description \code{sort_distributions} orders a given set of distribution assumptions
#' in order of intervals accuracy error in ascending or descending order.
#'
#' @param phi residual values used to compute the error interval.
#' @param errors set of real errors corresponding to the predictions of a particular
#' model.
#' @param dists character vector with the distribution assumptions to test. See also
#' 'Details'.
#' @param decreasing logical, indicating whether or not distributions should be 
#' ordered by decreasing accuracy error.
#' @param ... additional arguments to be passed to functions \code{error_interval}
#' and \code{acc_intervals}.
#'
#' @return Returns an object of class \code{c("df_intervals", "data.frame")} with
#' information of the distribution assumptions ordered by accuracy error.
#'
#' @details Allowed distribution assumptions are:
#' \itemize{
#'  \item{"n": }{Zero-mu Gaussian}
#'  \item{"nm": }{General Gaussian}
#'  \item{"l": }{Zero-mu Laplace}
#'  \item{"lm": }{General Laplace}
#'  \item{"b": }{Beta}
#'  \item{"w": }{Weibull}
#' }
#'
#' @seealso \link{df_intervals} \link{error_interval} \link{acc_intervals} \link{order}
#' @export
#'
#' @examples
#' sort_distributions(rnorm(10000),rnorm(10000),dists=c("n","b"))
#' sort_distributions(rnorm(10000),rnorm(10000),decreasing=TRUE)
sort_distributions <- function(phi,errors,dists=c("n","nm","l","lm","w","b")
                               ,decreasing=FALSE,...){
  distributions<-dists;
  errs<-vector();
  for (distribution in distributions){
    interv<-error_interval(phi, dist=distribution, ...);
    errs<-c(errs,acc_intervals(interv,errors,...)$err);
  }
  res<-df_intervals(distributions,errs);
  res<-res[order(res$error,decreasing=decreasing),];
  rownames(res)<-seq_along(res$error);
  return(res);
}
