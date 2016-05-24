#' Calibration sample data
#'
#' This includes the calibration sample data for the insurgency forecasting example in Montgomery, Hollenbach and Ward (2012). It provides the predictions for the three models included in the Ensemble model, as well as the true values of the dependent variable for insurgency in 29 Asian countries. The calibration sample ranges from January 2008 to December 2009. 
#'
#'
#' The variables included in the dataset are:
#' \itemize{
#' \item\code{LMER} The calibration sample predictions of the LMER model from the insurgency prediction example in Montgomery et. al. (2012). The LMER model is a generalized linear mixed effects model using the logistic link function. It includes two random effects terms and several other covariates.
#' \item\code{SAE} The calibration sample prediction of the SAE model from the insurgency prediction example in Montgomery et. al. (2012). This is a model developed as part of the ICEWS project and was designed by \emph{Strategic Analysis Enterprises}. It is a simple generalized linear model with 27 independent variables. 
#' \item\code{GLM} The calibration sample prediction of the GLM model from the insurgency prediction example in Montgomery et. al. (2012). This is a crude logistic model with only four independent variables. 
#' \item\code{Insurgency} The true values of the dependent variable in the calibration sample from the insurgency prediction example in Montgomery et. al. (2012). This is a binary variable indicating the actual ocurrence of insurgency for each observation in the calibration sample.
#'} 
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging. \emph{Political Analysis}. \emph{Political Analysis}. \bold{20}: 271-291.
#' @name calibrationSample
#' @docType data
NULL

#' Test sample data
#'
#' This includes the test sample data for the insurgency forecasting example in Montgomery, Hollenbach and Ward (2012). It provides the predictions for the three models included in the Ensemble model, as well as the true values of the dependent variable for insurgency in 29 Asian countries. The test sample ranges ranges from January 2010 to December 2010.  
#'
#' The variables included in the dataset are:
#' \itemize{
#' \item\code{LMER} The test sample predictions of the LMER model from the insurgency prediction example in Montgomery et. al. (2012). The LMER model is a generalized linear mixed effects model using the logistic link function. It includes two random effects terms and several other covariates.
#' \item\code{SAE} The test sample prediction of the SAE model from the insurgency prediction example in Montgomery et. al. (2012). This is a model developed as part of the ICEWS project and was designed by \emph{Strategic Analysis Enterprises}. It is a simple generalized linear model with 27 independent variables. 
#' \item\code{GLM} The test sample prediction of the GLM model from the insurgency prediction example in Montgomery et. al. (2012). This is a crude logistic model with only four independent variables. 
#' \item\code{Insurgency} The true values of the dependent variable in the test sample from the insurgency prediction example in Montgomery et. al. (2012). This is a binary variable indicating the actual ocurrence of insurgency for each observation in the test sample.
#'} 
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging.  \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @name testSample
#' @docType data
NULL

#' Sample data Presidential Election
#'
#' This includes the data for the presidential election forecasting example in Montgomery, Hollenbach and Ward (2012). The data ranges from 1952 to 2008 and includes predictions for the six different component models included in the Ensemble model. Users may split the sample into calibration and test sample. 
#'
#' The variables included in the dataset are:
#' \itemize{
#' \item\code{Campbell} Predictions of Campbell's ``Trial-Heat and Economy Model'' (Campbell 2008).
#' \item\code{Abramowitz} Predictions of Abramowitz's ``Time for Change Model'' (Abramowitz 2008).
#' \item\code{Hibbs} Predictions for the ``Bread and Peace Model'' created by Douglas Hibbs (2008).
#' \item\code{Fair} Forecasts from Fair's presidential vote share model (2010).  
#' \item\code{Lewis-Beck/Tien} Predictions from the ``Jobs Model Forecast''	by Michael Lewis-Beck and Charles Tien (2008).   
#' \item\code{EWT2C2} Predictions from the model in Column 2 in Table 2 by Erickson and Wlezien (2008). 
#' \item\code{Actual} The true values of the dependent variable, i.e. the incumbent-party voteshare in each presidential election in the sample.
#'} 
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging.  \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @references Campbell, James E. 2008. The trial-heat forecast of the 2008 presidential vote: Performance and value considerations in an open-seat election.  \emph{PS: Political Science & Politics} \bold{41}:697-701.
#'
#' @references Hibbs, Douglas A. 2000. Bread and peace voting in U.S presidential elections.  \emph{Public Choice} \bold{104}:149-180.
#'
#' @references Fair, Ray C. 2010. Presidential and Congressional vote-share equations: November 2010 update. Working Paper. Yale University.
#'
#' @references Lewis-Beck, Michael S. and Charles Tien. 2008. The job of president and the jobs model forecast: Obama for '08?  \emph{PS: Political Science & Politics} \bold{41}:687-690.
#'
#' @references Erikson, Robert S. and Christopher Wlezien. 2008. Leading economic indicators, the polls, and the presidential vote.  \emph{PS: Political Science & Politics} \bold{41}:703-707.

#' @name presidentialForecast
#' @docType data
NULL



#' EBMAforecast
#'
#' The EBMAforecast package (currently under development) allows users to increase the accuracy of forecasting models by pooling multiple component forecasts to generate ensemble forecasts. It includes functions to fit an ensemble Bayesian model averaging (EBMA) model using in-sample predictions, generate ensemble out-of-sample predictions, and create useful data visualizations.  Currently, the package can only handle dichotomous outcomes or those with normally distributed errors, although additional models will be added to the package in the coming months. Missing observation are allowed in the calibration set, but models with many predictions missing are penalized. 
#'
#' @name EBMAforecast
#' @docType package
#' @author  Michael D. Ward <\email{michael.d.ward@@duke.edu}> and Jacob M. Montgomery <\email{jacob.montgomery@@wustl.edu}> and Florian M. Hollenbach <\email{florian.hollenbach@@tamu.edu}>
#'
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2015). Calibrating ensemble forecasting models with sparse data in the social sciences.   \emph{International Journal of Forecasting}. In Press.
#' @references Montgomery, Jacob M., Florian M. Hollenbach and Michael D. Ward. (2012). Improving Predictions Using Ensemble Bayesian Model Averaging.   \emph{Political Analysis}. \bold{20}: 271-291.
#'
#' @references Raftery, A. E., T. Gneiting, F. Balabdaoui and M. Polakowski. (2005). Using Bayesian Model Averaging to calibrate forecast ensembles. \emph{Monthly Weather Review}. \bold{133}:1155--1174.
#' @references Sloughter, J. M., A. E. Raftery, T. Gneiting and C. Fraley. (2007). Probabilistic quantitative precipitation forecasting using Bayesian model averaging. \emph{Monthly Weather Review}. \bold{135}:3209--3220.
#' @references Fraley, C., A. E. Raftery, T. Gneiting. (2010). Calibrating Multi-Model Forecast Ensembles with Exchangeable and Missing Members using Bayesian Model Averaging. \emph{Monthly Weather Review}. \bold{138}:190--202.
#' @references Sloughter, J. M., T. Gneiting and A. E. Raftery. (2010). Probabilistic wind speed forecasting using ensembles and Bayesian model averaging. \emph{Journal of the American Statistical Association}. \bold{105}:25--35.
#' @references Fraley, C., A. E. Raftery, and T. Gneiting. (2010). Calibrating multimodel forecast ensembles with exchangeable and missing members using Bayesian model averaging. \emph{Monthly Weather Review}. \bold{138}:190--202.
#' @examples \dontrun{demo(EBMAforecast)
#' demo(presForecast)
#'}
NULL


