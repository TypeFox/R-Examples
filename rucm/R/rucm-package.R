#'rucm: Functions to model and predict a time series using Unobserved Components Model
#'
#'Package rucm contains functions to model and predict a time series using Unobserved Components Model (UCM) (Harvey (1989)) which decomposes the series into its salient components of trend, seasons, cycles, and regression effects due to predictors.
#'
#'Unobserved Components Models (UCMs) are special cases of more general and powerful tool in time series called State Space Models having an observation equation, which relates the dependent series to an unobserved state vector, and a state equation describing the evolution of the state vector over time. For a detailed discussion on State Space Models refer Harvey (1989) or Helske (2014). 
#'
#'
#'@references Harvey A. (1989). \emph{Forecasting, structural time series models and the Kalman filter}. Cambridge New York: Cambridge University Press
#'
#'Helske J (2014). \strong{KFAS}: \emph{Kalman filter and Smoothers for Exponential Family State Space Models}. R package version 1.0.4-1, URL \url{http://CRAN.R-project.org/package=KFAS}.
#'
#'SAS Institute Inc (2010). \emph{SAS/ETS 9.22 User's Guide}. SAS Institute Inc., Cary, NC. URL \url{http://support.sas.com/documentation/cdl/en/etsug/60372/PDF/default/etsug.pdf}.
#'
#'Selukar R (2011). "State Space Modeling Using SAS". \emph{Journal of Statistical Software}, \strong{41}(12), 1-13. URL \url{http://www.jstatsoft.org/v41/i12/}.
#'
#'Petris G, Petrone S (2011). "State Space Models in R". \emph{Journal of Statistical Software}, \strong{41}(4), 1-25. URL \url{http://www.jstatsoft.org/v41/i04/}.
#'
#'@examples
#'modelNile <- ucm(Nile~0, data = Nile, 
#'irregular = TRUE, level = TRUE, slope = TRUE)
#'
#'modelNile #Print the model
#'
#'#Return smoothed level values
#'modelNile$s.level 
#'
#'#Fixing the level variance to an absolute value
#'modelNile.fix <- ucm(Nile~0, data = Nile, 
#'irregular = TRUE, level = TRUE, level.var = 500, 
#'slope = TRUE) 
#'
#'#Predicting future values of the time series
#'predict(modelNile.fix, n.ahead = 12) 
#'
#'
#'@docType package
#'@name rucm
#'@rdname rucm-package
#'@aliases rucm
#'
NULL
