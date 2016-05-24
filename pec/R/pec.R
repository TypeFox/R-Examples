#' Prediction error curves
#' 
#' Evaluating the performance of risk prediction models in survival analysis.
#' The Brier score is a weighted average of the squared distances between the
#' observed survival status and the predicted survival probability of a model.
#' Roughly the weights correspond to the probabilities of not being censored.
#' The weights can be estimated depend on covariates. Prediction error curves
#' are obtained when the Brier score is followed over time.  Cross-validation
#' based on bootstrap resampling or bootstrap subsampling can be applied to
#' assess and compare the predictive power of various regression modelling
#' strategies on the same set of data.
#' 
#' Missing data in the response or in the input matrix cause a failure.
#' 
#' The status of the continuous response variable at cutpoints (\code{times}),
#' ie status=1 if the response value exceeds the cutpoint and status=0
#' otherwise, is compared to predicted event status probabilities which are
#' provided by the prediction models on the basis of covariates.  The
#' comparison is done with the Brier score: the quadratic difference between
#' 0-1 response status and predicted probability.
#' 
#' There are two different sources for bias when estimating prediction error in
#' right censored survival problems: censoring and high flexibility of the
#' prediction model. The first is controlled by inverse probability of
#' censoring weighting, the second can be controlled by special Monte Carlo
#' simulation. In each step, the resampling procedures reevaluate the
#' prediction model.  Technically this is done by replacing the argument
#' \code{object$call$data} by the current subset or bootstrap sample of the
#' full data.
#' 
#' For each prediction model there must be a \code{predictSurvProb} method: for
#' example, to assess a prediction model which evaluates to a \code{myclass}
#' object one defines a function called \code{predictSurvProb.myclass} with
#' arguments \code{object,newdata,cutpoints,...}
#' 
#' Such a function takes the object and
#' derives a matrix with predicted event status probabilities for each subject
#' in newdata (rows) and each cutpoint (column) of the response variable that
#' defines an event status.
#' 
#' Currently, \code{predictSurvProb} methods are available for the following
#' R-objects: \describe{ \item{}{\code{matrix}} \item{}{\code{aalen},
#' \code{cox.aalen} from \code{library(timereg)}} \item{}{\code{mfp} from
#' \code{library(mfp)}} \item{}{\code{phnnet}, \code{survnnet} from
#' \code{library(survnnet)}} \item{}{\code{rpart} (from \code{library(rpart)})}
#' \item{}{\code{coxph}, \code{survfit} from \code{library(survival)}}
#' \item{}{\code{cph}, \code{psm} from \code{library(rms)}}
#' \item{}{\code{prodlim} from \code{library(prodlim)}} \item{}{\code{glm} from
#' \code{library(stats)}} }
#'
#' @aliases pec
#' @param object A named list of prediction models, where allowed entries are
#' (1) R-objects for which a \link{predictSurvProb} method exists (see
#' details), (2) a \code{call} that evaluates to such an R-object (see
#' examples), (3) a matrix with predicted probabilities having as many rows as
#' \code{data} and as many columns as \code{times}. For cross-validation all
#' objects in this list must include their \code{call}.
#' @param formula A survival formula as obtained either
#' with \code{prodlim::Hist} or \code{survival::Surv}.
#' The left hand side is used to find the status response variable in \code{data}. For right censored data, the right
#' hand side of the formula is used to specify conditional censoring models.
#' For example, set \code{Surv(time,status)~x1+x2} and \code{cens.model="cox"}.
#' Then the weights are based on a Cox regression model for the censoring times
#' with predictors x1 and x2.  Note that the usual coding is assumed:
#' \code{status=0} for censored times and that each variable name that appears
#' in \code{formula} must be the column name in \code{data}. If there are no
#' covariates, i.e. \code{formula=Surv(time,status)~1} the \code{cens.model} is
#' coerced to \code{"marginal"} and the Kaplan-Meier estimator for the
#' censoring times is used to calculate the weights.  If \code{formula} is
#' \code{missing}, try to extract a formula from the first element in object.
#' @param data A data frame in which to validate the prediction models and to
#' fit the censoring model.  If \code{data} is missing, try to extract a data
#' set from the first element in object.
#' @param traindata A data frame in which the models are trained. This argument
#' is used only in the absence of crossvalidation, in which case it is
#' passed to the predictHandler function predictSurvProb 
#' @param times A vector of time points. At each time point the prediction
#' error curves are estimated. If \code{exact==TRUE} the \code{times} are
#' merged with all the unique values of the response variable.  If \code{times}
#' is missing and \code{exact==TRUE} all the unique values of the response
#' variable are used.  If missing and \code{exact==FALSE} use a equidistant
#' grid of values between \code{start} and \code{maxtime}.  The distance is
#' determined by \code{exactness}.
#' @param cause For competing risks, the event of interest. Defaults to the
#' first state of the response, which is obtained by evaluating the left hand
#' side of \code{formula} in \code{data}.
#' @param start Minimal time for estimating the prediction error curves.  If
#' missing and \code{formula} defines a \code{Surv} or \code{Hist} object then
#' \code{start} defaults to \code{0}, otherwise to the smallest observed value
#' of the response variable. \code{start} is ignored if \code{times} are given.
#' @param maxtime Maximal time for estimating the prediction error curves. If
#' missing the largest value of the response variable is used.
#' @param exact Logical. If \code{TRUE} estimate the prediction error curves at
#' all the unique values of the response variable. If \code{times} are given
#' and \code{exact=TRUE} then the \code{times} are merged with the unique
#' values of the response variable.
#' @param exactness An integer that determines how many equidistant gridpoints
#' are used between \code{start} and \code{maxtime}.  The default is 100.
#' @param fillChar Symbol used to fill-in places where the values of the
#' prediction error curves are not available. The default is \code{NA}.
#' @param cens.model Method for estimating inverse probability of censoring
#' weigths:
#' 
#' \code{cox}: A semi-parametric Cox proportional hazard model is fitted to the
#' censoring times
#' 
#' \code{marginal}: The Kaplan-Meier estimator for the censoring times
#' 
#' \code{nonpar}: Nonparametric extension of the Kaplan-Meier for the censoring
#' times using symmetric nearest neighborhoods -- available for arbitrary many
#' strata variables on the right hand side of argument \code{formula} but at
#' most one continuous variable. See the documentation of the functions
#' \code{prodlim} and \code{neighborhood} from the prodlim package.
#' 
#' \code{aalen}: The nonparametric Aalen additive model fitted to the censoring
#' times. Requires the \code{timereg} package.
#' @param ipcw.refit If \code{TRUE} the inverse probability of censoring
#' weigths are estimated separately in each training set during
#' cross-validation.
#' @param ipcw.args List of arguments passed to function specified by argument \code{cens.model}.
#' @param splitMethod SplitMethod for estimating the prediction error curves.
#' 
#' \code{none/noPlan}: Assess the models in the same data where they are
#' fitted.  \code{boot}: DEPRECIATED.
#' 
#' \code{cvK}: K-fold cross-validation, i.e. \code{cv10} for 10-fold
#' cross-validation.  After splitting the data in K subsets, the prediction
#' models (ie those specified in \code{object}) are evaluated on the data
#' omitting the Kth subset (training step). The prediction error is estimated
#' with the Kth subset (validation step).
#' 
#' The random splitting is repeated \code{B} times and the estimated prediction
#' error curves are obtained by averaging.
#' 
#' \code{BootCv}: Bootstrap cross validation. The prediction models are trained
#' on \code{B} bootstrap samples, that are either drawn with replacement of the
#' same size as the original data or without replacement from \code{data} of
#' the size \code{M}.  The models are assessed in the observations that are NOT
#' in the bootstrap sample.
#' 
#' \code{Boot632}: Linear combination of AppErr and BootCvErr using the
#' constant weight .632.
#' 
#' \code{Boot632plus}: Linear combination of AppErr and BootCv using weights
#' dependent on how the models perform in permuted data.
#' 
#' \code{loocv}: Leave one out cross-validation.
#' 
#' \code{NoInf}: Assess the models in permuted data.
#' @param B Number of bootstrap samples. The default depends on argument
#' \code{splitMethod}.  When \code{splitMethod} in
#' c("BootCv","Boot632","Boot632plus") the default is 100. For
#' \code{splitMethod="cvK"} \code{B} is the number of cross-validation cycles,
#' and -- default is 1.  For \code{splitMethod="none"} \code{B} is the number
#' of bootstrap simulations e.g. to obtain bootstrap confidence limits --
#' default is 0.
#' @param M The size of the bootstrap samples for resampling without
#' replacement. Ignored for resampling with replacement.
#' @param reference Logical. If \code{TRUE} add the marginal Kaplan-Meier
#' prediction model as a reference to the list of models.
#' @param model.args List of extra arguments that can be passed to the
#' \code{predictSurvProb} methods. The list must have an entry for each entry
#' in \code{object}.
#' @param model.parms Experimental.  List of with exactly one entry for each
#' entry in \code{object}.  Each entry names parts of the value of the fitted
#' models that should be extracted and added to the value.
#' @param keep.index Logical. If \code{FALSE} remove the bootstrap or
#' cross-validation index from the output list which otherwise is included in
#' the splitMethod part of the output list.
#' @param keep.matrix Logical. If \code{TRUE} add all \code{B} prediction error
#' curves from bootstrapping or cross-validation to the output.
#' @param keep.models Logical. If \code{TRUE} keep the models in object. Since
#' fitted models can be large objects the default is \code{FALSE}.
#' @param keep.residuals Logical. If \code{TRUE} keep the patient individual
#' residuals at \code{testTimes}.
#' @param keep.pvalues For \code{multiSplitTest}. If \code{TRUE} keep the
#' pvalues from the single splits.
#' @param noinf.permute If \code{TRUE} the noinformation error is approximated
#' using permutation.
#' @param multiSplitTest If \code{TRUE} the test proposed by van de Wiel et al.
#' (2009) is applied. Requires subsampling bootstrap cross-validation, i.e.
#' that \code{splitMethod} equals \code{bootcv} and that \code{M} is specified.
#' @param testIBS A range of time points for testing differences between models
#' in the integrated Brier scores.
#' @param testTimes A vector of time points for testing differences between
#' models in the time-point specific Brier scores.
#' @param confInt Experimental.
#' @param confLevel Experimental.
#' @param verbose if \code{TRUE} report details of the progress, e.g. count the
#' steps in cross-validation.
#' @param savePath Place in your file system (i.e., a directory on your
#' computer) where training models fitted during cross-validation are saved. If
#' \code{missing} training models are not saved.
#' @param slaveseed Vector of seeds, as long as \code{B}, to be given to the
#' slaves in parallel computing.
#' @param na.action Passed immediately to model.frame. Defaults to na.fail. If
#' set otherwise most prediction models will not work.
#' @param ... Not used.
#' @return A \code{pec} object. See also the help pages of the corresponding
#' \code{print}, \code{summary}, and \code{plot} methods.  The object includes
#' the following components: \item{PredErr}{ The estimated prediction error
#' according to the \code{splitMethod}. A matrix where each column represents
#' the estimated prediction error of a fit at the time points in time.  }
#' \item{AppErr}{ The training error or apparent error obtained when the
#' model(s) are evaluated in the same data where they were trained. Only if
#' \code{splitMethod} is one of "NoInf", "cvK", "BootCv", "Boot632" or
#' "Boot632plus".  } \item{BootCvErr}{ The prediction error when the model(s)
#' are trained in the bootstrap sample and evaluated in the data that are not
#' in the bootstrap sample.  Only if \code{splitMethod} is one of "Boot632" or
#' "Boot632plus". When \code{splitMethod="BootCv"} then the \code{BootCvErr} is
#' stored in the component \code{PredErr}.  } \item{NoInfErr}{ The prediction
#' error when the model(s) are evaluated in the permuted data.  Only if
#' \code{splitMethod} is one of "BootCv", "Boot632", or "Boot632plus".  For
#' \code{splitMethod="NoInf"} the \code{NoInfErr} is stored in the component
#' \code{PredErr}.  } \item{weight}{ The weight used to linear combine the
#' \code{AppErr} and the \code{BootCvErr} Only if \code{splitMethod} is one of
#' "Boot632", or "Boot632plus".  } \item{overfit}{ Estimated \code{overfit} of
#' the model(s).  See Efron \& Tibshirani (1997, Journal of the American
#' Statistical Association) and Gerds \& Schumacher (2007, Biometrics).  Only
#' if \code{splitMethod} is one of "Boot632", or "Boot632plus".  }
#' \item{call}{The call that produced the object} \item{time}{The time points
#' at which the prediction error curves change.} \item{ipcw.fit}{The fitted
#' censoring model that was used for re-weighting the Brier score residuals.
#' See Gerds \& Schumacher (2006, Biometrical Journal)} \item{n.risk}{The
#' number of subjects at risk for all time points.} \item{models}{The
#' prediction models fitted in their own data.} \item{cens.model}{The censoring
#' models.} \item{maxtime}{The latest timepoint where the prediction error
#' curves are estimated.} \item{start}{The earliest timepoint where the
#' prediction error curves are estimated.} \item{exact}{\code{TRUE} if the
#' prediction error curves are estimated at all unique values of the response
#' in the full data.} \item{splitMethod}{The splitMethod used for estimation of
#' the overfitting bias.}
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{plot.pec}}, \code{\link{summary.pec}},
#' \code{\link{R2}}, \code{\link{crps}}
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' 
#' E. Graf et al.  (1999), Assessment and comparison of prognostic
#' classification schemes for survival data. Statistics in Medicine, vol 18,
#' pp= 2529--2545.
#' 
#' Efron, Tibshirani (1997) Journal of the American Statistical Association 92,
#' 548--560 Improvement On Cross-Validation: The .632+ Bootstrap Method.
#' 
#' Gerds, Schumacher (2006), Consistent estimation of the expected Brier score
#' in general survival models with right-censored event times. Biometrical
#' Journal, vol 48, 1029--1040.
#' 
#' Thomas A. Gerds, Martin Schumacher (2007) Efron-Type Measures of Prediction
#' Error for Survival Analysis Biometrics, 63(4), 1283--1287
#' doi:10.1111/j.1541-0420.2007.00832.x
#' 
#' Martin Schumacher, Harald Binder, and Thomas Gerds. Assessment of survival
#' prediction models based on microarray data. Bioinformatics, 23(14):1768-74,
#' 2007.
#' 
#' Mark A. van de Wiel, Johannes Berkhof, and Wessel N. van Wieringen Testing
#' the prediction error difference between 2 predictors Biostatistics (2009)
#' 10(3): 550-560 doi:10.1093/biostatistics/kxp011
#' @keywords survival
#' @examples
#' 
#' # simulate an artificial data frame
#' # with survival response and two predictors
#' 
#' set.seed(130971)
#' library(prodlim)
#' library(survival)
#' dat <- SimSurv(100)
#' 
#' # fit some candidate Cox models and compute the Kaplan-Meier estimate 
#' 
#' Models <- list("Cox.X1"=coxph(Surv(time,status)~X1,data=dat,y=TRUE),
#'               "Cox.X2"=coxph(Surv(time,status)~X2,data=dat,y=TRUE),
#'               "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat,y=TRUE))
#' 
#' # compute the apparent prediction error 
#' PredError <- pec(object=Models,
#'                   formula=Surv(time,status)~X1+X2,
#'                   data=dat,
#'                   exact=TRUE,
#'                   cens.model="marginal",
#'                   splitMethod="none",
#'                   B=0,
#'                   verbose=TRUE)
#' 
#' print(PredError,times=seq(5,30,5))
#' summary(PredError)
#' plot(PredError,xlim=c(0,30))
#' 
#' # Comparison of Weibull model and Cox model
#' library(survival)
#' library(rms)
#' library(pec)
#' data(pbc)
#' pbc <- pbc[sample(1:NROW(pbc),size=100),]
#' f1 <- psm(Surv(time,status!=0)~edema+log(bili)+age+sex+albumin,data=pbc)
#' f2 <- coxph(Surv(time,status!=0)~edema+log(bili)+age+sex+albumin,data=pbc)
#' f3 <- cph(Surv(time,status!=0)~edema+log(bili)+age+sex+albumin,data=pbc,surv=TRUE)
#' brier <- pec(list("Weibull"=f1,"CoxPH"=f2,"CPH"=f3),data=pbc,formula=Surv(time,status!=0)~1)
#' print(brier)
#' plot(brier)
#' 
#' # compute the .632+ estimate of the generalization error
#' set.seed(130971)
#' library(prodlim)
#' library(survival)
#' dat <- SimSurv(100)
#' set.seed(17100)
#' PredError.632plus <- pec(object=Models,
#'                   formula=Surv(time,status)~X1+X2,
#'                   data=dat,
#'                   exact=TRUE,
#'                   cens.model="marginal",
#'                   splitMethod="Boot632plus",
#'                   B=3,
#'                   verbose=TRUE)
#' 
#' print(PredError.632plus,times=seq(4,12,4))
#' summary(PredError.632plus)
#' plot(PredError.632plus,xlim=c(0,30))
#' # do the same again but now in parallel
#' \dontrun{set.seed(17100)
#' library(doMC)
#' registerDoMC()
#' PredError.632plus <- pec(object=Models,
#'                   formula=Surv(time,status)~X1+X2,
#'                   data=dat,
#'                   exact=TRUE,
#'                   cens.model="marginal",
#'                   splitMethod="Boot632plus",
#'                   B=3,
#'                   verbose=TRUE)
#' }
#' # assessing parametric survival models in learn/validation setting
#' learndat <- SimSurv(50)
#' testdat <- SimSurv(30)
#' library(rms)
#' f1 <- psm(Surv(time,status)~X1+X2,data=learndat)
#' f2 <- psm(Surv(time,status)~X1,data=learndat)
#' pf <- pec(list(f1,f2),formula=Surv(time,status)~1,data=testdat,maxtime=200)
#' plot(pf)
#' summary(pf)
#' 
#' # ---------------- competing risks -----------------
#' 
#' library(survival)
#' library(riskRegression)
#' library(cmprsk)
#' library(pec)
#' cdat <- SimCompRisk(100)
#' data(cdat)
#' f1  <- CSC(Hist(time,event)~X1+X2,cause=2,data=cdat)
#' f1a  <- CSC(Hist(time,event)~X1+X2,cause=2,data=cdat,survtype="surv")
#' f2  <- CSC(Hist(time,event)~X1,data=cdat,cause=2)
#' require(cmprsk)
#' ## predict.crr <- cmprsk:::predict.crr
#' f3  <- FGR(Hist(time,event)~X1+X2,cause=2,data=cdat)
#' f4  <- FGR(Hist(time,event)~X1+X2,cause=2,data=cdat)
#' p1 <- pec(list(f1,f1a,f2,f3,f4),formula=Hist(time,event)~1,data=cdat,cause=2)
#' 
#' @export 
# {{{ header pec.list
pec <- function(object,
                formula,
                data,
                traindata,
                times,
                cause,
                ## time points
                start,
                maxtime,
                exact=TRUE,
                exactness=100,
                fillChar=NA,
                ## censoring weighting 
                cens.model="cox",
                ipcw.refit=FALSE,
                ipcw.args=NULL,
                ## data splitting
                splitMethod="none",
                B,
                M,
                ## misc parameters
                reference=TRUE,
                model.args=NULL,
                model.parms=NULL,
                keep.index=FALSE,
                keep.matrix=FALSE,
                keep.models=FALSE,
                keep.residuals=FALSE,
                keep.pvalues=FALSE,
                noinf.permute=FALSE,
                multiSplitTest=FALSE,
                testIBS,
                testTimes,
                confInt=FALSE,
                confLevel=0.95,
                verbose=TRUE,
                savePath=NULL,
                slaveseed=NULL,
                na.action=na.fail,
                ...)
{
    # }}}
    # {{{ checking integrity some arguments

  theCall=match.call()
  if (match("replan",names(theCall),nomatch=FALSE))
    stop("The argument name 'replan' has been replaced by 'splitMethod'.")
  
  if (!missing(testIBS) && (!(is.logical(testIBS) || (length(testIBS)==2 && is.numeric(testIBS)))))
    stop("Argument testIBS can be TRUE/FALSE or a vector of two numeric values.")
  if (missing(testIBS)) testIBS <- FALSE
  if (keep.residuals && missing(testTimes))
    stop("To keep.residuals please specify testTimes.")
  if (missing(splitMethod) && multiSplitTest==TRUE){
    stop("Need data splitting to compute van de Wiel's test")
  }
  if (missing(M) && multiSplitTest) M <- NA

  # }}}
  # {{{ check and convert object
  if (class(object)[1]!="list") {
      object <- list(object)
  }
  # }}}
  # {{{ formula
  if (missing(formula)){
      if (length(grep("~",as.character(object[[1]]$call$formula)))==0){
          stop(paste("Argument formula is missing and first model has no usable formula:",as.character(object[[1]]$call$formula)))
      } else{
          ftry <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
          if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
              stop("Argument formula is missing and first model has no usable formula.")
          else if (verbose)
              warning("Formula missing. Using formula from first model")
      }
  }
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
      stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
  }
  else{
      if (!(formula.names[2] %in% c("Surv","Hist")))
          survp <- FALSE
      else
          survp <- TRUE
  }
  # }}}
  # {{{ data
  if (missing(data)){
      data <- eval(object[[1]]$call$data)
      if (match("data.frame",class(data),nomatch=0)==0)
          stop("Argument data is missing.")
      else
          if (verbose)
              warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  # }}}
  # {{{ censoring model
  
  cens.model <- match.arg(cens.model,c("cox","marginal","nonpar","aalen","none","rfsrc"))
  
  # }}}
  # {{{ response

  histformula <- formula
  if (histformula[[2]][[1]]==as.name("Surv")){
      histformula[[2]][[1]] <- as.name("Hist")
  }
  ## m <- model.frame(histformula,data,na.action=na.fail)
  m <- model.frame(histformula,data,na.action=na.action)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=0)!=0){
      attr(response,"model") <- "survival"
      attr(response,"cens.type") <- "rightCensored"
      model.type <- "survival"
  }
  model.type <- attr(response,"model")
  if (model.type=="competing.risks"){
      predictHandlerFun <- "predictEventProb"
      if (missing(cause))
          cause <- attr(response,"state")[1]
  }
  else{
      if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
      if (survp==TRUE && NCOL(response)!=2) stop("Survival response must have two columns: time and status.")
      predictHandlerFun <- "predictSurvProb"
  }
  
  # }}}
  # {{{ prediction models
  if (reference==TRUE) {
      ProdLimform <- as.formula(update(formula,".~NULL"))
      ## ProdLimfit <- do.call(prodlim::prodlim,list(formula=ProdLimform,data=data))
      ProdLimfit <- prodlim::prodlim(formula=ProdLimform,data=data)
      ProdLimfit$call$data <- as.character(substitute(data))
      ProdLimfit$call$formula=ProdLimform
      ProdLimfit$formula <- as.formula(ProdLimfit$formula)
      ## print(environment(ProdLimfit$formula))
      ## if (model.type=="competing.risks")
      ## object <- c(list(Reference=ProdLimfit),object)
      ## else
      ## browser()
      object <- c(list("Reference"=ProdLimfit),object)
  }
  if (is.null(names(object))){
      names(object) <- sapply(object,function(o)class(o)[1])
      names(object) <- make.names(names(object),unique=TRUE)
  }
  else{ # fix missing names
      if (any(names(object)=="")){
          names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
          names(object) <- make.names(names(object),unique=TRUE)
      }else{
          # leave names as they were given
      }
  }
  ## names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object)
  # }}}  
  # {{{ sort the data 

  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
    if (predictHandlerFun=="predictEventProb"){
      event <- prodlim::getEvent(response,mode="character")
      event <- event[neworder]
    }
    response <- response[neworder,,drop=FALSE]
    Y <- response[,"time"]
    status <- response[,"status"]
  }
  else{
    cens.model <- "none"
    neworder <- order(response)
    Y <- response[neworder]
    status <- rep(1,length(Y))
  }
  ## for competing risks find the cause of interest.
  if (predictHandlerFun=="predictEventProb"){
    availableCauses <- unique(event)
    if (!match(cause, availableCauses,nomatch=FALSE))
      stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    event <- event==cause
  }
  ##   else{
  ##     event <- NULL
  ##   }
  data <- data[neworder,]
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)
  # }}}
  # {{{ splitMethod

  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  B <- splitMethod$B
  ResampleIndex <- splitMethod$index
  k <- splitMethod$k
  do.resample <- !(is.null(ResampleIndex))
  if (keep.matrix==TRUE & !do.resample){
    warning("Argument keep.matrix set to FALSE, since no resampling/crossvalidation is requested.")
    keep.matrix <- FALSE
  }

  # }}}      
  # {{{ find maxtime, start, and jumptimes in the range of the response 
if (missing(maxtime) || is.null(maxtime))
    maxtime <- unique.Y[NU]
  if (missing(start))
    if (survp==TRUE)
      start <- 0  ## survival times are positive
    else
      start <- min(unique.Y) 
  if (missing(times)){
    if (exact==TRUE)
      times <- unique(c(start,unique.Y))
    else
      times <- seq(start,maxtime,(maxtime - start)/exactness)
  }
  else{
    if (exact==TRUE) 
      times <- sort(c(start,unique(times),unique.Y))
    else
      times <- sort(unique(c(start,times)))
  }
  times <- times[times<=maxtime]
  NT <-  length(times)

# }}}
# {{{ IPCW (all equal to 1 without censoring) 

  if((cens.model %in% c("aalen","cox","nonpar"))){
      if (all(as.numeric(status)==1) || sum(status)==N){
          if (verbose)
              message("No censored observations: cens.model coerced to \"none\".")
          cens.model <- "none"
      }
      if ((cens.model!="nonpar") && length(attr(terms(formula),"factors"))==0){
          if (verbose==TRUE)
              message("No covariates  specified: Kaplan-Meier for censoring times used for weighting.")
          cens.model <- "marginal"}
  }
  if (predictHandlerFun=="predictEventProb"){
      iFormula <- as.formula(paste("Surv(itime,istatus)","~",as.character(formula)[[3]]))
      iData <- data
      iData$itime <- response[,"time"]
      iData$istatus <- response[,"status"]
      if (ipcw.refit==TRUE)
          stop("pec: internal refitting of censoring distribution not (not yet) supported for competing risks")
      ipcw.call <- NULL
      ipcw <- ipcw(formula=iFormula,
                   data=iData,
                   method=cens.model,
                   args=ipcw.args,
                   times=times,
                   subjectTimes=Y,
                   subjectTimesLag=1)
      ipcw$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
  }
  else{
      if (ipcw.refit==TRUE && splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632"))
          ipcw.call <- list(formula=formula,data=NULL,method=cens.model,times=times,subjectTimes=NULL,subjectTimesLag=1)
      else
          ipcw.call <- NULL
      ipcw <- ipcw(formula=formula,
                   data=data,
                   method=cens.model,
                   args=ipcw.args,
                   times=times,
                   subjectTimes=Y,
                   subjectTimesLag=1)
      ipcw$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
  }
  ## force ipc weights not to exaggerate
  ## weights should not be greater than 1/(sample size)
  ## if (ipcw$dim==1){
  ## ipcw$IPCW.times <- apply(ipcw$IPCW.times,1,function(x)pmax(x,1/N))
  ## }
  ## else{
  ## ipcw$IPCW.times <- pmax(ipcw$IPCW.times,1/N)
  ## }
  ## ipcw$IPCW.subjectTimes <- pmax(ipcw$IPCW.subjectTimes,1/N)
  ## browser()
  
  #  wt <- ipcw$IPCW.times
  #  wt.obs <- ipcw$IPCW.subjectTimes
  #  if (NCOL(wt)>1) {stopifnot(length(wt)==(N*NT))}  else{stopifnot(length(wt)==NT)}

  # }}}
# {{{ checking the models for compatibility with resampling
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  # }}}
# {{{ ---------------------------Apparent error---------------------------

  AppErr <- lapply(1:NF,function(f){
      ## message(f)
      fit <- object[[f]]
      extraArgs <- model.args[[f]]
      if (predictHandlerFun=="predictEventProb"){ # competing risks
          pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,cause=cause),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,,drop=FALSE]
      .C("pecCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
    }
    else{  # survival
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,,drop=FALSE]
      .C("pec",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
    }
  })
  names(AppErr) <- names(object)

  ## se.Apperr <- lapply(1:NF,function(f){
  ## ## message(f)
  ## fit <- object[[f]]
  ## extraArgs <- model.args[[f]]
  ## if (predictHandlerFun=="predictEventProb"){ # competing risks
  ## pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,cause=cause),extraArgs))
  ## if (class(object[[f]])[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,,drop=FALSE]
  ## Paulo(as.double(Y),
  ## as.double(status),
  ## as.double(event),
  ## as.double(times),
  ## as.double(pred),
  ## as.double(ipcw$IPCW.times),
  ## as.double(ipcw$IPCW.subjectTimes))
  ## }
  ## else{  # survival
  ## pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times),extraArgs))
  ## if (class(object[[f]])[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,,drop=FALSE]
  ## Paulo(as.double(Y),
  ## as.double(status),
  ## as.double(times),
  ## as.double(pred),
  ## as.double(ipcw$IPCW.times),
  ## as.double(ipcw$IPCW.subjectTimes))
  ## }})


  # }}}
# {{{------------------------No information error------------------------

  if (splitMethod$internal.name %in% c("Boot632plus")){
    if (verbose==TRUE){
      message("Computing noinformation error using all permutations")
    }
    if (noinf.permute==FALSE){
      NoInfErr <- lapply(1:NF,function(f){
        fit <- object[[f]]
        extraArgs <- model.args[[f]]
        pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times),extraArgs))
        extraArgs <- model.args[[f]]
        if (predictHandlerFun=="predictEventProb")
          .C("pec_noinfCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
        else
          .C("pec_noinf",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
      })
      names(NoInfErr) <- names(object)
    }else{
      if (verbose==TRUE){
        message("Noinformation error simulation loop (B=",B,")")
      }
      ## FIXME: need to parallelize noinf
      NoInfErrList <- lapply(1:B,function(b){
        if (verbose==TRUE){
          internalTalk(b,B,sign=".")
        }
        responseNames <- colnames(response)
        noinf.b <- data[sample(1:NROW(data),replace=FALSE),-match(responseNames,names(data))]
        noinf.b[,responseNames] <- response
        ipcw.b <- ipcw(formula=formula,data=noinf.b,method=cens.model,args=ipcw.args,times=times,subjectTimes=Y,subjectTimesLag=1)
        noinfPredErr <- lapply(1:NF,function(f){
          fit.b <- internalReevalFit(object=object[[f]],data=noinf.b,step=b,silent=FALSE,verbose=verbose)
          ## fit.b$call <- object[[f]]$call
          extraArgs <- model.args[[f]]

          pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times),extraArgs))
          if (predictHandlerFun=="predictEventProb"){
            pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times,cause=cause),extraArgs))
            .C("pecCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred.b),as.double(ipcw.b$IPCW.times),as.double(ipcw.b$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
          }
          else{
            pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times),extraArgs))
            .C("pec",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred.b),as.double(ipcw.b$IPCW.times),as.double(ipcw.b$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
          }
        })
        noinfPredErr
      })
      NoInfErrMat <- lapply(1:NF,function(f){
        do.call("rbind",lapply(NoInfErrList,function(x){
          x[[f]]
        }))})
      NoInfErr <- lapply(NoInfErrMat,colMeans)
      names(NoInfErr) <- names(object)
    }
  }

  # }}}
# {{{--------------k-fold and leave-one-out CrossValidation-----------------------

  if (splitMethod$internal.name %in% c("crossval","loocv")){
      kCV <- kFoldCrossValidation(object=object,data=data,Y=Y,status=status,event=event,times=times,cause=cause,ipcw=ipcw,splitMethod=splitMethod,giveToModel=model.args,predictHandlerFun=predictHandlerFun,keep=keep.matrix,verbose=verbose)
      CrossValErr <- kCV$CrossValErr
      if (keep.matrix && B>1)
          CrossValErrMat <- kCV$CrossValErrMat
  }

  # }}}
# {{{ ----------------------BootstrapCrossValidation----------------------

  if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
    if (verbose==TRUE){
      message("Split sample loop (B=",B,")")
    }
    if (missing(testTimes)){
      testTimes <- NULL
    }
    BootCv <- bootstrapCrossValidation(object=object,
                                       data=data,
                                       Y=Y,
                                       status=status,
                                       event=event,
                                       times=times,
                                       cause=cause,
                                       ipcw=ipcw,
                                       ipcw.refit=ipcw.refit,
                                       ipcw.call=ipcw.call,
                                       splitMethod=splitMethod,
                                       multiSplitTest=multiSplitTest,
                                       testIBS=testIBS,
                                       testTimes=testTimes,
                                       confInt=confInt,
                                       confLevel=confLevel,
                                       getFromModel=model.parms,
                                       giveToModel=model.args,
                                       predictHandlerFun=predictHandlerFun,
                                       keepMatrix=keep.matrix,
                                       keepResiduals=keep.residuals,
                                       verbose=verbose,
                                       savePath=savePath,
                                       slaveseed=slaveseed)
    BootstrapCrossValErr <- BootCv$BootstrapCrossValErr
    Residuals <- BootCv$Residuals
    names(BootstrapCrossValErr) <- names(object)
    if (multiSplitTest==TRUE){
      comparisons <- allComparisons(names(object))
      multiSplitTestResults <- list(testIBS=testIBS,B=B,M=M,N=N,testTimes=testTimes)
      multiSplitTestResults$Comparisons <- lapply(1:length(comparisons),function(cc){
        if (length(testTimes)>0){
          allPairwisePvaluesTimes <- do.call("rbind",lapply(BootCv$testedResid,function(b){
            b$pValue[[cc]]}))
          out <- list(pValueTimes=apply(allPairwisePvaluesTimes,2,median))
          if (keep.pvalues==TRUE){
            out$allPairwisePvaluesTimes <- allPairwisePvaluesTimes
          }
        }
        else out <- NULL
        if(length(testIBS)>0){
          allPairwisePvaluesIBS <- sapply(BootCv$testedResid,function(b){
            b$IBSpValue[[cc]]
          })
          out$pValueIBS <- median(allPairwisePvaluesIBS)
        }
        if (keep.pvalues==TRUE){
          out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS}
        out
      })
      names(multiSplitTestResults$Comparisons) <- names(comparisons)
      ## multiSplitTest$splitMethod <- splitMethod
      class(multiSplitTestResults) <- "multiSplitTest"
    }
    ## upperLimits <- lapply(BootCv$testedResid,function(x){x[,1:length(testTimes)]})
    ##     if (testIBS==TRUE){
    ##       wtestIBSpValues <- do.call("cbind",apply(BootCv$testedResid,function(x){x[,length(testTimes)+1]}))
    ##     }
    ## wtestIBSupper <- BootCv$testedResid$wtestIBSupper
    ##   }
    if (keep.matrix==TRUE){
      BootstrapCrossValErrMat <- BootCv$BootstrapCrossValErrMat
      names(BootstrapCrossValErr) <- names(object)
    }
  }

  # }}}
  # {{{ Bootstrap .632
  if (splitMethod$internal.name=="Boot632"){
      B632Err <- lapply(1:NF,function(f){
          .368 * AppErr[[f]] + .632 * BootstrapCrossValErr[[f]]
      })
      names(B632Err) <- names(object)
  }
  # }}}    
  # {{{ Bootstrap .632+

  if (splitMethod$internal.name=="Boot632plus"){
    B632plusErr <- lapply(1:NF,function(f){
      Err1 <- pmin(BootstrapCrossValErr[[f]],NoInfErr[[f]])
      overfit <- (Err1 - AppErr[[f]]) / (NoInfErr[[f]] - AppErr[[f]])
      overfit[!(Err1>AppErr[[f]])] <- 0
      w <- .632 / (1 - .368 * overfit)
      B632plusErr <- (1-w) * AppErr[[f]]  + w * Err1
      B632plusErr
      ## w[NoInfErr<=BootstrapCrossValErr] <- 1
      ## B632plus.error <- (1-w) * AppErr  + w * BootstrapCrossValErr
    })
    names(B632plusErr) <- names(object)
  }

  # }}}
# {{{ prepare output

out <- switch(splitMethod$internal.name,
              "noPlan"=list("AppErr"=AppErr),
              "Boot632plus"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"NoInfErr"=NoInfErr,"Boot632plusErr"=B632plusErr),
              "Boot632"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"Boot632Err"=B632Err),
              "BootCv"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr),
              "loocv"=list("AppErr"=AppErr,"loocvErr"=CrossValErr),
              "crossval"=list("AppErr"=AppErr,"crossvalErr"=CrossValErr),
              "noinf"=list("AppErr"=AppErr,"NoInfErr"=NoInfErr))
observed.maxtime <- sapply(out,function(x){
    ## lapply(x,function(y){times[length(y)-sum(is.na(y))-1]})
    lapply(x,function(y){times[length(y)-sum(is.na(y))]})
})
minmaxtime <- min(unlist(observed.maxtime))
if (multiSplitTest==TRUE){
    out <- c(out,list(multiSplitTest=multiSplitTestResults))
}
if (keep.residuals==TRUE){
    out <- c(out,list(Residuals=Residuals))
}
if (keep.matrix==TRUE && splitMethod$internal.name!="noPlan"){
    if (splitMethod$internal.name %in% c("crossval","loocv")){
        if (B>1)
            out <- c(out,list("CrossValErrMat"=CrossValErrMat))
    }
    else{
        if (splitMethod$internal.name!="noinf")
            out <- c(out,list("BootstrapCrossValErrMat"=BootstrapCrossValErrMat))
    }
}
if (!is.na(fillChar))
    out <- lapply(out,function(o){
        o[is.na(o)] <- fillChar
        o
    })
if (!is.null(model.parms))
    out <- c(out,list("ModelParameters"=BootCv$ModelParameters))
  
  
  if (!keep.index) splitMethod$index <- NULL
  n.risk <- N - prodlim::sindex(Y,times)
  # }}}
  # {{{ put out
  if(keep.models==TRUE){
      outmodels <- object
  } else{
        outmodels <- names(object)
        names(outmodels) <- names(object)
    }
  out <- c(out,
           list(call=theCall,
                response=model.response(m),
                time=times,
                ## ipcw.fit=as.character(ipcw$fit$call),
                n.risk=n.risk,
                models=outmodels,
                maxtime=maxtime,
                observed.maxtime=observed.maxtime,
                minmaxtime=minmaxtime,
                reference=reference,
                start=min(times),
                cens.model=cens.model,
                exact=exact,
                splitMethod=splitMethod))
  ##   if (verbose==TRUE && splitMethod$internal.name %in% c("BootCv","Boot632","Boot632plus","crossval","loocv")) cat("\n")
  class(out) <- "pec"
  out

  # }}}
}


