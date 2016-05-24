#' This function maximizes an arbitrary likelihood including generalized linear models and Cox partial likelihood.
#'
#' This function provides flexible model fitting. The main model specification is written in \code{modelexpr}. \code{modelexpr} consisits of two parts separated by \sQuote{~}. The distribution is specified in the first part, and the second part includes variable name which follows the specified distribution in the first part. Available distributional specifications are \sQuote{cox}(Cox partial likelihood), \sQuote{pois}(Poisson distribution), \sQuote{norm}(normal distribution), \sQuote{binom}(binomial distribution), \sQuote{nbinom}(negative binomial distribution), \sQuote{gamma}(gamma distribution) and \sQuote{weibull}(Weibull distribution). Options can be specified for some distribution specification after \sQuote{/} in the distribution specification part. One optional specification format is \sQuote{optionname=value}. Multiple options separated by \sQuote{,} can also be specified.
#'
#' For Cox regressions, time and status variable must be specified in parentheses like \code{cox(time, status)}. Some options are also available for Cox regressions, and \sQuote{efron}, \sQuote{breslow}, \sQuote{discrete} and \sQuote{average} is available for tie handling method. \code{ties(discrete)} specification corresponds to \sQuote{exact} ties specification in \code{coxph} function, and \code{ties(average)} corresonds to \sQuote{exact} specification in SAS PHREG procecure. See references for further details. Strata option which specifies a variable indicating strata is also available in Cox regressions. Subset option which has same functinality as subset argument below is also available for Cox regressions and other distribution specifications. For other distribution specifications, parameters must be specified in parentheses. For poisson distribution, mean parameter must be specified as \code{pois(mean)}. Note that each parameter specificaiton can be a variable or R equation. For other distributions, \code{norm(mean, variance)}, \code{binom(size, probability)}, \code{nbinom(size, probability)}, \code{gamma(shape, scale)}, \code{weibull(shape, scale)}.
#'
#' When distributions are specified, additional R expressions can be specified in \code{preexpr} argument. R expressions are parsed to make variable list. Variables which exist in data.frame or the global environment must be vector, and the rest of variables are regarded as parameters. If you define variable \sQuote{mu} in \code{preexpr}, you can use \sQuote{mu} in \code{modelexpr} argument. Refer Poisson regression examples below.
#'
#' Time dependent covariate in Cox regression is supported experimentally. Time can be referred as \code{time_inner_} at event time. Refer time dependent covariate example below.
#' 
#' @title Model fitting function for epifit package
#' @param modelexpr a character string or string vector specifying the model. See \sQuote{Details}.
#' @param preexpr a character string or string vector specifying R expressions executed before \sQuote{modelexpr}. Multiple expressions separated by \sQuote{;} is allowed. See \sQuote{Details}.
#' @param timedepexpr a character string or string vector specifying R expressions executed during Cox regression at each time of event occurrence. Typical use is incorporating time-dependent variables. Not yet implemented completely.
#' @param nullparam a numeric vector specifying null model. The default is zero value is streched to a vector of same length as parameters, which indicate a null model where all parameter values are equal to zero.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument.
#' @param subset an expression indicating which subset of the rows of data should be used in the fit. All observations are included by default.
#' @param weight vector of case weights. If weights is a vector of integers, then the estimated coefficients are equivalent to estimating the model from data with the individual cases replicated as many times as indicated by weights. Only supported \sQuote{breslow} and \sQuote{efron} ties specification in Cox regression models and other likelihood specifications.
# @param random a character string specifying random effects. See \sQuote{Details}.
#' @param na.action a missing-data filter function. This is applied when data.frame is supplied as \sQuote{data} parameter. Default is \code{options()$na.action}.
#' @param opt a character string specifying the method for optimization. When sQuote{newrap} is specified, \code{nlm} function that uses Newton-type algorithm used for obtaining maximum likelihood estimate. For the rest of specifications (\sQuote{BFGS} for quasi-Newton method, \sQuote{CG} for conjugate gradients method, \sQuote{SANN} for simulated annealing, \sQuote{Nelder-Mead} for Nelder and Mead simplex method), \code{optim} is used. The default is \sQuote{newrap}.
#' @param tol1 a numeric value specifying \code{gtol} in nlm, \code{abstol} in optim.
#' @param tol2 a numeric value specifying \code{stol} in nlm, \code{reltol} in optim.
#' @param maxiter a integer value specifying the maximum number of iterations. Defaults is 200.
#' @param init a numeric vector or list specifying initial values for the parameters specified in the form of vector. To search for the best initial values for parameters, provide a numeric vector including candidate values for a parameter with the name of its parameter. A default value can be specified as a list component without name. For example, there are three parameters (beta0, beta1, beta2) in a model and specified as \code{init=list(0, beta1=c(0,1,2)}, initial parameter values tried is (0, 0, 0), (0, 1, 0), and (0, 2, 0), and the best initial value is used for estimation.
#' @param verbatim a integer value from 0 (minimum) to 2 (maximum) controlling the amount of information printed during calculation.
#' @param ... for the arguments used in the inner functions (currently not used).
#' @return a list containing the result of model fitting including parameter estimates, variance of parameter estimates, log likelihood and so on.
#' @useDynLib epifit Rf_select Rf_mumul
#' @references DeLong, D. M., Guirguis, G.H., and So, Y.C. (1994). Efficient computation of subset selection probabilities with application to Cox regression. \emph{Biometrika} \strong{81}, 607-611.
#' @references Gail, M. H., Lubin, J. H., and Rubinstein, L. V. (1981). Likelihood calculations for matched case-control studies and survival studies with tied death times. \emph{Biometrika} \strong{68}, 703-707.
#' @seealso
#' \code{\link{AIC.epifit}},
#' \code{\link{print.epifit}}
#' @examples
#' library(survival)
#' 
#' # Make sample data
#' set.seed(123)
#' nsub <- 20
#' follow <- 5
#' x <- rnorm(nsub)
#' rate <- exp(-0.5 + x)
#' etime <- rweibull(nsub, 1, 1/rate)
#' status <- as.integer(etime < follow)
#' time <- pmin(follow, etime)
#' dat <- data.frame(status, x, time)
#'
#' coxph(Surv(time, status)~x, data=dat)
#' modelexpr <- "cox(time,status)~exp(beta*x)"
#' epifit(modelexpr=modelexpr, data=dat)
#'
#' glm(status ~ x + offset(log(time)), family=poisson(), data=dat)
#' preexpr <- "mu <- time*exp(beta0 + beta1*x)"
#' modelexpr <- "pois(mu) ~ status"
#' epifit(modelexpr=modelexpr, preexpr=preexpr, data=dat)
#'
#' # The simplest test data set from coxph function
#' test1 <- list(time=c(4,3,1,1,2,2,3),
#'               status=c(1,1,1,0,1,1,0),
#'               x=c(0,2,1,1,1,0,0),
#'               sex=c(0,0,0,0,1,1,1))
#'
#' # Cox regressions with strata
#' coxph(Surv(time, status) ~ x + strata(sex), data=test1)
#' modelexpr <- "cox(time,status)/strata(sex)~exp(beta*x)"
#' epifit(modelexpr=modelexpr, data=test1)
#'
#' # Tie specification example in Cox regressions
#' coxph(Surv(time, status) ~ x + strata(sex), data=test1, ties="breslow")
#' modelexpr <- "cox(time,status)/strata=sex,ties=breslow~exp(beta*x)"
#' epifit(modelexpr=modelexpr, data=test1)
#' 
#' # Average partial likelihood
#' modelexpr <- "cox(time,status)/strata=sex,ties=average~exp(beta*x)"
#' epifit(modelexpr=modelexpr, data=test1)
#' 
#' # Conditional logistic regression for matched case-control studies
#' # hypothetical data
#' conlog <- data.frame(strata=c(1,1,2,2,3,3,4,4,5,5), outcome=c(1,0,1,0,1,0,1,0,1,0),
#'                      cov=c(1,3,2,1,5,2,4,2,2,2))
#' # Make dummy survival time so that all the cases in a matched set have the same survival
#' # time value, and the corresponding controls are censored at later times
#' conlog <- cbind(conlog, dummy=(2 - conlog$outcome))
#' coxph(Surv(dummy, outcome)~cov + strata(strata), ties="exact", data=conlog)
#' modelexpr <- "cox(dummy,outcome)/ties=discrete,strata=strata~exp(beta*cov)"
#' epifit(modelexpr=modelexpr, data=conlog)
#'
#' 
#' # Joint model example (for demonstrating technical specifications)
#' nsub <- 1000
#' follow <- 10
#' x <- rnorm(nsub)
#' dose <- rweibull(nsub, 0.5, 1/(2.84)^2)
#' rate <- exp(-1 + x)*(1 + 0.5*dose)
#'
#' # Generate survival data
#' etime <- rweibull(nsub, 1, 1/rate)
#' status <- as.integer(etime < follow)
#' time <- pmin(follow, etime)
#' dat2 <- data.frame(event=status, py=time, x, dose, model=1)
#'
#' # Generate person-year table (baseline is different)
#' py <- runif(nsub)*follow
#' rate2 <- exp(-0.5 + 0.5*x)*(1 + 0.5*dose)
#' event <- sapply(rate2*py, function(x){rpois(1, x)})
#' dat3 <- cbind(pytable(event, py, cbind(x,dose)), model=2)
#' dat4 <- rbind(dat2, dat3)
#'
#' modelexpr <- c("cox(py,event)/subset=(model==1)~exp(beta0*x)*(1 + beta*dose)",
#'              "pois(py*exp(beta1 + beta2*x)*(1 + beta*dose))/subset=(model==2) ~ event")
#' epifit(modelexpr, data=dat4)
#'
#' # Time dependent covariate example
#' id <- 1:8
#' group <- c(0, 0, 0, 0, 1, 1, 1, 1)
#' time <- c(4, 5, 7, 9, 6, 10, 11, 12)
#' event <- c(1, 1, 0, 1, 1, 1, 1, 0)
#' dat5 <- data.frame(id=id, group=group, time=time, event=event)
#' modelexpr <- "cox(time, event) ~ exp(beta1*group + beta2*t_g)"
#' # t_g is time-dependent variable created by using event time time_inner_ (created automatically)
#' timedepexpr <- "t_g <- time_inner_ * group"
#' epifit(modelexpr=modelexpr, timedepexpr=timedepexpr, data=dat5)
#' @export
epifit <- function(modelexpr=NULL, preexpr="", timedepexpr="", nullparam=NULL,
                   data=NULL, subset=NULL, weight=NULL, na.action=NULL,
                   opt=c("newrap", "BFGS", "CG", "Nelder-Mead", "SANN"),
                   tol1=1e-8, tol2=1e-8, maxiter=200, init=NULL, verbatim=0, ...){

  ## argument processing
  args <- list(...)
  Call <- match.call(expand.dots=TRUE)
  opt <- match.arg(opt)

  if(is.null(modelexpr))
    stop("modelexpr cannot be omitted")

  # scalar variable
  num_nmodel <- length(modelexpr) # Num of models in modelexpr
  num_npreexpr <- numeric(num_nmodel) # Num of sentence in preexpr
  num_ntimedepexpr <- numeric(num_nmodel) # number of sentence in preexpr
  num_nsubject <- numeric(1) # number of data
  num_nullparam <- numeric(1) # number of parameters to estimate in null model
  
  ## variable list for all model
  vec_parameter <- NULL
  vec_variable <- NULL
  vec_varparam <- NULL # variance parameter list for null model

  ## null parameter position (to be estimated)
  vec_nullpos <- NULL

  ## model name vector
  vec_modelname <- character(num_nmodel)

  ## inner variable name (weight, time1, time2, status, time)
  vec_innername=c("weight_inner_", "time1_inner_", "time2_inner_",
    "status_inner_", "time_inner_")

  ## parameter init value candidate
  lst_initparam <- list()
  
  ## parsed symbols
  psd_preexpr <- list()
  psd_timedepexpr <- list()

  ## model information per model
  lst_option <- list() # option list
  lst_psdparams <- list() # parsed 1st, 2nd 3rd element of modelexpr *******

  ## variable information per model
  lst_parameter <- list() # parameter list per model (do not exist in the data)
  lst_variable <- list() # variable list per model (exist in the data)
  lst_assigned <- list() # assigned variable list per model (non exist, but created)
  lst_dependent <- list() # dependent variable list (almost same as lst_asssigned)
  
  if(is.null(na.action)){
    na.action <- options()[["na.action"]]
  }

  if(na.action=="na.fail"){
    
    if(is.null(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.fail(data), parent=.GlobalEnv)
    }

  } else if(na.action=="na.omit"){

    if(is.null(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.omit(data), parent=.GlobalEnv)
    }
    
  } else if(na.action=="na.exclude"){

    if(is.null(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.exclude(data), parent=.GlobalEnv)
    }
    
  } else if(na.action=="na.pass"){

    if(is.null(data)){
      data <- .GlobalEnv
    } else {
      data <- list2env(na.pass(data), parent=.GlobalEnv)
    }

  } else {
      stop("invalid na.action specification")
  }

  ## one model
  if(num_nmodel == 1){

    tryCatch(
      psd_preexpr[[1]] <- parse(text=preexpr),
      error=function(e){stop("Error in parsing preexpr argument")}
      )
    tryCatch(
      psd_timedepexpr[[1]] <- parse(text=timedepexpr),
      error=function(e){stop("Error in parsing timedepexpr argument")}
      )

    if(max(length(psd_preexpr), length(timedepexpr)) > num_nmodel)
      warning("One of the number of preexpr or timedepexpr is larger than that of modelexpr")

    ## more than one model
  } else {

    ## preexpr
    if(length(preexpr) != num_nmodel){
      if(length(preexpr) != 1){
        stop("number of preexpr is not equal to that of models")
      }
      # in case of one expression for all models
      tryCatch(
        for(i in 1:num_nmodel){psd_preexpr[[i]] <- parse(text=preexpr)},
        error=function(e){stop("Error in parsing preexpr argument")}
        )
    } else {
      for(i in 1:num_nmodel)
        tryCatch(
          psd_preexpr[[i]] <- parse(text=preexpr[i]),
          error=function(e){stop("Error in parsing preexpr argument")}
          )
    }

    ## timedepexpr
    if(length(timedepexpr) != num_nmodel){
      if(length(timedepexpr) != 1){
        stop("number of timedepexpr is not equal to that of models")
      }
      # in case of one expression for all models
      tryCatch(
        for(i in 1:num_nmodel){psd_timedepexpr[[i]] <- parse(text=timedepexpr)},
        error=function(e){stop("Error in parsing timedepexpr argument")}
        )     
    } else {
      for(i in 1:num_nmodel)
        tryCatch(
          psd_timedepexpr[[i]] <- parse(text=timedepexpr[i]),
          error=function(e){stop("Error in parsing timedepexpr argument")}
          )
    }
  }

  for(i in 1:num_nmodel){
    
    lst_option[[i]] <- list() # option list
    lst_psdparams[[i]] <- list() # parsed 1st, 2nd 3rd element of modelexpr *******

    ## per model
    lst_parameter[[i]] <- list() # parameter list per model (do not exist in the data)
    lst_variable[[i]] <- list() # variable list per model (exist in the data)
    lst_assigned[[i]] <- list() # assigned variable list per model (non exist, but created)
    lst_dependent[[i]] <- list() # dependent variable list (almost same as lst_asssigned)
  
    ## left term is response variable specification
    if(length(grep(";", modelexpr[i])) != 0)
      stop("\";\" cannot be included in modelexpr")

    tmpstr <- strsplit(modelexpr[i], "~")[[1]]
    
    ## left side of "modelexpr" with option specifications
    tmp_leftmodel <- trim(tmpstr[1])
    
    ## right term of "modelexpr" should be evaluated
    tmp_rightmodel <- trim(tmpstr[2])
    
    lst_option[[i]] <- GetOptions(tmp_leftmodel)
    tmp_leftmodel<- strsplit(tmp_leftmodel, "/")[[1]][1] # remove option specification
    
    ## parse vec_response
    tryCatch(
      psd_leftmodel <- parse(text=tmp_leftmodel)[[1]],
      error=function(e){stop("Error in parsing modelexpr[", i, "] argument")}
      )    

    tryCatch(
      psd_rightmodel <- parse(text=tmp_rightmodel)[[1]],
      error=function(e){stop("Error in parsing modelexpr[", i, "] argument")}
      )    

    ## model name
    vec_modelname[i] <- as.character(psd_leftmodel[[1]])

    if(length(tmpstr) != 2){
      stop("\"modelexpr\" must include response variable as left term separated by \"~\"")
    }
    
    ## obtain number of sentences in each expression
    num_npreexpr[i] <- length(psd_preexpr[[i]])
    num_ntimedepexpr[i] <-  length(psd_timedepexpr[[i]])
    
    ## check number of parameters 
    if(length(psd_leftmodel) < 2){
      stop("parameter should be specified in modelexpr")
    } else if(length(psd_leftmodel) == 2){
      if(!vec_modelname[i] %in% c("pois", "general")){
        stop("number of parameters specified in \"modelexpr\" is incorrect")
      }
      lst_psdparams[[i]] <- list(psd_rightmodel, psd_leftmodel[[2]])
    } else if(length(psd_leftmodel) == 3){
      if(!vec_modelname[i] %in% c("cox", "norm", "binom", "nbinom", "gamma", "weibull")){
        stop("number of parameters specified in \"modelexpr\" is incorrect")
      }
      lst_psdparams[[i]] <- list(psd_rightmodel, psd_leftmodel[[2]], psd_leftmodel[[3]])
    } else if(length(psd_leftmodel) == 4){
      if(!vec_modelname[i] %in% c("cox")){
        stop("number of parameters specified in \"modelexpr\" is incorrect")
      }
      lst_psdparams[[i]] <- list(psd_rightmodel, psd_leftmodel[[2]], psd_leftmodel[[3]], psd_leftmodel[[4]])
    } else {
        stop("number of parameters specified in \"modelexpr\" is incorrect")
    }
    ## lst_psdparams, lst_option are finished

    for(j in 1:length(lst_psdparams[[i]])){
      tmp <- ListVariable(lst_psdparams[[i]][[j]])
      lst_assigned[[i]] <- c(lst_assigned[[i]], tmp[[1]])
      tmp <- ClassifyParameter(tmp[[2]], data, tmp[[1]])
      lst_parameter[[i]] <- c(lst_parameter[[i]], tmp[[1]])
      lst_variable[[i]] <- c(lst_variable[[i]], tmp[[2]])
      if(j==3 && vec_modelname[i]=="norm"){
        ## make variance parameter list
        vec_varparam <- c(vec_varparam, unlist(tmp[[1]]))
      }
    }

    if(num_npreexpr[i] > 0){
      tmp <- ListVariable(psd_preexpr[[i]])
      for(j in 1:length(tmp)){
        lst_assigned[[i]] <- c(lst_assigned[[i]], tmp[[i]][[1]])
        tmp2 <- ClassifyParameter(tmp[[i]][[2]], data, tmp[[i]][[1]])
        lst_parameter[[i]] <- c(lst_parameter[[i]], tmp2[[1]])
        lst_variable[[i]] <- c(lst_variable[[i]], tmp2[[2]])
      }
    }

    if(num_ntimedepexpr[i] > 0){
      tmp <- ListVariable(psd_timedepexpr[[i]])
      for(j in 1:length(tmp)){
        lst_assigned[[i]] <- c(lst_assigned[[i]], tmp[[i]][[1]])
        tmp2 <- ClassifyParameter(tmp[[i]][[2]], data, tmp[[i]][[1]])
        lst_parameter[[i]] <- c(lst_parameter[[i]], tmp2[[1]])
        lst_variable[[i]] <- c(lst_variable[[i]], tmp2[[2]])
      }
    }

    ## obtain dependent variable
    lst_dependent[[i]] <- LimitVarlist(unlist(lst_parameter[[i]]), unlist(lst_assigned[[i]]))
  }

  ## summarize parameters and variables 
  vec_parameter <- unique(RemoveVariableName(unlist(lst_parameter),
                                             c(unlist(lst_dependent), vec_innername)))
  vec_variable <- unique(unlist(lst_variable))

  if(length(vec_parameter) < 1){
    stop("there are no parameters to estimate")
  }
  
  ## obtain number of subjects from the first variable
  num_nsubject <- length(get(vec_variable[1], envir=data))

  if(!is.null(weight)){
    if(is.character(weight)){
      weight <- get(weight, envir=data)
    }
  }

  ## Get intial value from init argument
  if(!is.null(init)){
    if(is.numeric(init)){
      init <- rep(init, length.out=length(vec_parameter))
      
    } else if(is.list(init)){
      initparamname <- names(init)
      if(is.null(initparamname)){
        stop("List specification in init argument must be provided with name")
      } else {
        pos <- match("", initparamname, nomatch=0)
        defaultval <- init[[pos]][1] # only use first element
        if(!is.numeric(defaultval))
          stop("List element of init must be numeric value or vector")
        for(i in 1:length(vec_parameter)){
          pos <- match(vec_parameter[i], initparamname, nomatch=0)
          if(pos > 0){
            lst_initparam[[i]] <- init[[pos]]
            if(!is.numeric(defaultval))
              stop("List element of init must be numeric value or vector")
          } else {
            lst_initparam[[i]] <- defaultval
          }
        }
        names(lst_initparam) <- vec_parameter
      }
      init <- do.call(makeInitVector, lst_initparam)

    } else {
      warning("Argument \"init\" must be numeric vector or list, and replaced with default values")
      init <- rep(0, length(vec_parameter))
    }
  } else {
    init <- rep(0, length(vec_parameter))
    # set 1 to variance parameter
    if(length(vec_varparam) > 0)
      init[GetParamPosition(vec_varparam, vec_parameter)] <- 1
  }

  envs <- list(NULL)
  envstrata <- list(NULL)
  
  ## make environments envs[[model, strata, random]]
  for(i in 1:num_nmodel){

    ## subset processing (commented out for debug)
    if(!is.null(subset)){
      if(is.character(subset))
        tryCatch(
          index <- rep(TRUE, num_nsubject) & as.logical(eval(parse(text=subset), envir=data)),
          error=function(e){stop("Error in executing subset argument")}
          )      
      else if(is.numeric(subset)){
        index <- (subset != 0)
      } else if(is.logical(subset)) {
        index <- subset
      }
      if(length(subset) != num_nsubject){
        warning("Length of subset argument is different from data length")
        subset <- subset & rep(TRUE, num_nsubject)
      }
    } else {
      index <- rep(TRUE, num_nsubject)
    }

    ## Cox regression
    if(vec_modelname[i]=="cox"){
      
      if(!is.null(lst_option[[i]][["subset"]])){

        tryCatch(
          index <- index & as.logical(eval(parse(text=lst_option[[i]][["subset"]]), envir=data)),
          error=function(e){stop("Error in executing subset option")}
          )
      }

      if(length(lst_psdparams[[i]]) == 4){

        tryCatch(
          time1 <- eval(parse(text=lst_psdparams[[i]][[1]]), envir=data),
          error=function(e){stop("Error in evaluating argument in Cox regression")}
          )

        tryCatch(
          time2 <- eval(parse(text=lst_psdparams[[i]][[2]]), envir=data),
          error=function(e){stop("Error in evaluating argument in Cox regression")}
          )
        
        tryCatch(
          status <- eval(parse(text=lst_psdparams[[i]][[3]]), envir=data),
          error=function(e){stop("Error in evaluating argument in Cox regression")}
          )

      } else if(length(lst_psdparams[[i]]) == 3){

        time1 <- rep(0, num_nsubject)

        tryCatch(
          time2 <- eval(parse(text=lst_psdparams[[i]][[2]]), envir=data),
          error=function(e){stop("Error in evaluating argument in Cox regression")}
          )
        
        tryCatch(
          status <- eval(parse(text=lst_psdparams[[i]][[3]]), envir=data),
          error=function(e){stop("Error in evaluating argument in Cox regression")}
          )
        
      } else {
        stop("unsupported model specification for Cox regression")
      }

      ## status range check
      if(min(status[index]) < 0 || max(status[index]) > 1){
        stop("range of status variable is not betweenn 0 and 1")
      }

      if(length(time1) != length(time2) || length(time2) != length(status) ||
         length(status) != num_nsubject)
        stop("Length of time or status variable in Cox regression is different")

      ## clear dim information
      dim(time1) <- NULL
      dim(time2) <- NULL
      dim(status) <- NULL
                
      ## Make environments (strata/random variable)
      if(!is.null(lst_option[[i]][["strata"]])){ # with strata

        tryCatch(
          strata <- eval(parse(text=lst_option[[i]][["strata"]]), envir=data)[index],
          error=function(e){stop("Error in evaluating strata option in Cox regression")}
          )        

        stratalist <- unique(strata[order(strata)])
        
        for(j in 1:length(stratalist)){
          
          index2 <- index & (strata==stratalist[j]) # strata restriction within model
          
          { # BLOCK
            envstrata[[j]] <- new.env()
            orderedtime <- order(time2[index2])

            for(k in 1:length(vec_variable)){
              assign(vec_variable[k], get(vec_variable[k], envir=data)[index2][orderedtime], envir=envstrata[[j]])
            }

            if(!is.null(weight))
              assign(vec_innername[1], weight[index2][orderedtime], envir=envstrata[[j]])
  
            assign(vec_innername[2], time1[index2][orderedtime], envir=envstrata[[j]])
            assign(vec_innername[3], time2[index2][orderedtime], envir=envstrata[[j]])
            assign(vec_innername[4], status[index2][orderedtime], envir=envstrata[[j]])
            
          } # BLOCK end
        }
        
      } else { # without strata

        { # BLOCK
          envstrata[[1]] <- new.env()
          orderedtime <- order(time2[index])
          for(k in 1:length(vec_variable)){
            assign(vec_variable[k], get(vec_variable[k], envir=data)[index][orderedtime], envir=envstrata[[1]])
          }

          if(!is.null(weight))
            assign(vec_innername[1], weight[index][orderedtime], envir=envstrata[[1]])
          
          assign(vec_innername[2], time1[index][orderedtime], envir=envstrata[[1]])
          assign(vec_innername[3], time2[index][orderedtime], envir=envstrata[[1]])
          assign(vec_innername[4], status[index][orderedtime], envir=envstrata[[1]])
          
        } # BLOCK end
      }
      
      envs[[i]] <- envstrata
      
    ## other likelihood models      
    } else {
      
      if(!is.null(lst_option[[i]][["strata"]])){
        warning("strata option is ignored because supported only in Cox regression model")
      }
      
      if(!is.null(lst_option[[i]][["subset"]])){
        index <- index & eval(parse(text=lst_option[[i]][["subset"]]), envir=data)
      }

      { # No random BLOCK
        
        envs[[i]] <- list(new.env())

        CopyVariableBetweenEnvironment(vec_variable, envs[[i]][[1]], data, index, TRUE) 
        
        if(!is.null(weight))
          assign(vec_innername[1], weight[index], envir=envs[[i]][[1]])
        
      } # BLOCK end
    }
  }

  ## Calculate best initial value
  if(is.list(init)){

    WrapperLoglikelihood <- function(...){
      tryCatch(return(LogLikelihood(...)),
               error=function(e){return(NA)}
               )
    }

    PrintVerbatim(0, verbatim, length(init),"\nNumber of initial parameter combinations searched for is ")

    initresult <- sapply(init, WrapperLoglikelihood,
                         num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
                         lst_option=lst_option, lst_psdparams=lst_psdparams,
                         psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
                         vec_parameter=vec_parameter, vec_variable=vec_variable,
                         vec_innername=vec_innername
                         )
    initresult <- cbind(1:length(init), initresult)
    initresult <- initresult[order(initresult[,2]),]
    init <- init[[initresult[1,1]]]
  }

  ## Obtain null model
  if(is.null(nullparam)){
    nullparam <- rep(0, length.out=length(vec_parameter))
    pos <- unlist(GetParamPosition(vec_varparam, vec_parameter))
    if(length(pos) > 0){
      nullparam[pos] <- NA
    }
  } else {
    if(length(nullparam) != length(vec_parameter))
      stop("Length of nullparam must be equal to the number of parameters")
  }
  
  vec_nullpos <- is.na(nullparam)
  num_nullparam <- sum(is.na(nullparam))

  result <- NULL
  
  names(init) <- vec_parameter
  
  PrintVerbatim(0, verbatim, init, "\nChoice of initial value of parameters are:\n")

  if(num_nullparam == 0){
  
    names(nullparam) <- vec_parameter
    PrintVerbatim(1, verbatim, nullparam, "\nNull parameter values are:\n")

    tryCatch(
      nulllik <- LogLikelihood(init=nullparam,
                               num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
                               lst_option=lst_option, lst_psdparams=lst_psdparams,
                               psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
                               vec_parameter=vec_parameter, vec_variable=vec_variable,
                               vec_innername=vec_innername),
      error=function(e){print(e);
                        stop("Error in evaluating loglikelihood function at null value")}
      )
  } else {

    WrapperLogLikelihood <- function(init, nullparam, vec_nullpos, ...){
      nullparam[vec_nullpos] <- init
      LogLikelihood(init=nullparam, ...)
    }

    if(opt=="newrap"){
      
      tryCatch(
        ans <- nlm(f=WrapperLogLikelihood, nullparam=nullparam, vec_nullpos=vec_nullpos,
                   p=init[vec_nullpos], iterlim=maxiter, hessian=TRUE, gradtol=tol1,
                   steptol=tol2, print.level=verbatim,
                   num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
                   lst_option=lst_option, lst_psdparams=lst_psdparams,
                   psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
                   vec_parameter=vec_parameter, vec_variable=vec_variable,
                   vec_innername=vec_innername),
        error=function(e){stop("Error in calculating loglikelihood for nullmodel")}
        )
      
      nulllik <- ans$minimum
      
    } else {

      tryCatch(
        ans <- optim(par=init[vec_nullpos], fn=WrapperLogLikelihood, nullparam=nullparam,
                     vec_nullpos=vec_nullpos, method=opt, hessian=TRUE,
                     control=list(reltol=tol1, maxit=maxiter, trace=as.integer(verbatim==2),
                       abstol=tol2),
                     num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
                     lst_option=lst_option, lst_psdparams=lst_psdparams,
                     psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
                     vec_parameter=vec_parameter, vec_variable=vec_variable,
                     vec_innername=vec_innername),
        error=function(e){stop("Error in calculating loglikelihood for nullmodel")}
        )
      
      nulllik <- ans$value
      
    }
  }
  
  result <- list(n=num_nsubject, call=Call, modelformula=modelexpr,
                 parameters=vec_parameter, variables=vec_variable)
  
  class(result) <- "epifit"
  attr(result, "df") <- length(vec_parameter)
  
  if(opt=="newrap"){
    
    ans <- nlm(f=LogLikelihood, p=init, iterlim=maxiter, hessian=TRUE, gradtol=tol1,
               steptol=tol2, print.level=verbatim,
               num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
               lst_option=lst_option, lst_psdparams=lst_psdparams,
               psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
               vec_parameter=vec_parameter, vec_variable=vec_variable,
               vec_innername=vec_innername)
    
    result <- MakeResultFromNlm(result, ans, nulllik)
    
  } else {
    
    ans <- optim(par=init, fn=LogLikelihood, method=opt, hessian=TRUE,
                 control=list(reltol=tol1, maxit=maxiter, trace=as.integer(verbatim==2),
                              abstol=tol2),
                 num_nmodel=num_nmodel, vec_modelname=vec_modelname, envs=envs,
                 lst_option=lst_option, lst_psdparams=lst_psdparams,
                 psd_preexpr=psd_preexpr, psd_timedepexpr=psd_timedepexpr,
                 vec_parameter=vec_parameter, vec_variable=vec_variable,
                 vec_innername=vec_innername)
    
    result <- MakeResultFromOptim(result, ans, nulllik)
    
  }

  return(result)
}
