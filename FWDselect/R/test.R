#'Bootstrap based test for covariate selection
#'@description Function that applies a bootstrap based test for covariate
#'  selection. It helps to determine the number of variables to be included in
#'  the model.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param method A character string specifying which regression method is used,
#'  i.e., linear models (\code{"lm"}), generalized additive models.
#'@param family A description of the error distribution and link function to be
#'  used in the model: (\code{"gaussian"}), (\code{"binomial"}) or
#'  (\code{"poisson"}).
#'@param nboot Number of bootstrap repeats.
#'@param speedup A logical value. If  \code{TRUE} (default), the testing procedure
#' is  computationally efficient since it considers one more variable to fit
#' the alternative model than the number of variables used to fit the null.
#' If \code{FALSE}, the fit of the alternative model is based on considering
#' the best subset of variables of size greater than \code{q}, the one that minimizes an
#' information criterion. The size of this subset must be given by the user
#' filling the argument \code{qmin}.
#' @param qmin By default \code{NULL}. If \code{speedup} is \code{FALSE},
#' \code{qmin} is an integer number selected by the user. To help you select
#' this argument, it is recommended to visualize the graphical output of the
#' \code{plot} function and choose the number \code{q} which minimizes the curve.
#'@param unique A logical value. By default \code{FALSE}. If  \code{TRUE},
#' the test is performed only for one null hypothesis, given by the argument  \code{q}.
#'@param q By default \code{NULL}. If  \code{unique} is \code{TRUE}, \code{q}
#' is the size of the subset of variables to be tested.
#'@param bootseed Seed to be used in the bootstrap procedure.
#'@param cluster A logical value. If  \code{TRUE} (default), the testing
#'  procedure is  parallelized.
#'@param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores to be used
#' is equal to the number of cores of the machine - 1.
#'@details In a regression framework, let \eqn{X_1, X_2, \ldots, X_p},  a set of
#'  \eqn{p} initial variables and \eqn{Y} the response variable, we propose a
#'  procedure to test the null hypothesis of \eqn{q} significant variables  in
#'  the model --\eqn{q}  effects not equal to zero-- versus the alternative in
#'  which the model contains more than \eqn{q} variables. Based on the general
#'  model \deqn{Y=m(\textbf{X})+\varepsilon  \quad {\rm{where}} \quad
#'  m(\textbf{X})= m_{1}(X_{1})+m_{2}(X_{2})+\ldots+m_{p}(X_{p})} the following
#'  strategy is considered: for a subset of size \eqn{q}, considerations will be
#'  given to a test for the null hypothesis \deqn{H_{0} (q): \sum_{j=1}^p
#'  I_{\{m_j \ne 0\}} \le q} vs. the general hypothesis \deqn{H_{1} :
#'  \sum_{j=1}^p I_{\{m_j \ne 0\}} > q}
#'@return
#' A list with two objects. The first one is a table containing
#'\item{Hypothesis}{Number of the null hypothesis tested}
#'\item{Statistic}{Value of the T statistic}
#'\item{pvalue}{pvalue obtained in the testing procedure}
#'\item{Decision}{Result of the test for a significance level of 0.05}
#'
#'The second argument \code{nvar} indicates the number of variables that
#' have to be included in the model.
#'@references Sestelo, M., Villanueva, N. M. and Roca-Pardinas, J. (2013).
#'  FWDselect: an R package for selecting variables in regression models.
#'  Discussion Papers in Statistics and Operation Research, University of Vigo, 13/01.
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@note The detailed expression of the formulas are described in HTML help
#'  \url{http://cran.r-project.org/web/packages/FWDselect/FWDselect.pdf}
#'@seealso \code{\link{selection}}
#'@examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' test(x, y, method = "lm", cluster = FALSE, nboot = 5)
#'
#' ## for speedup = FALSE
#' # obj2 = qselection(x, y, qvector = c(1:9), method = "lm",
#' # cluster = FALSE)
#' # plot(obj2) # we choose q = 7 for the argument qmin
#' # test(x, y, method = "lm", cluster = FALSE, nboot = 5,
#' # speedup = FALSE, qmin = 7)
#'
#'@importFrom parallel detectCores
#'@importFrom parallel clusterExport
#'@importFrom parallel parCapply
#'@importFrom parallel stopCluster
#'@importFrom stats as.formula
#'@importFrom stats lm
#'@importFrom stats glm
#'@importFrom stats predict
#'@importFrom stats rbinom
#'@importFrom stats rpois
#'@export

test <- function(x, y, method = "lm", family = "gaussian", nboot = 50,
                 speedup = TRUE, qmin = NULL, unique = FALSE, q = NULL,
                 bootseed = NULL, cluster = TRUE, ncores = NULL){

  # Statistics T
  Tvalue <- function(xy, qT = qh0, optionT = method,
                     speed = speedup, prevars = NULL) {
    x = xy[, 2:ncol(xy)]
    y = xy[, 1]
    var_res = NULL
    nvar = ncol(x)
    x = as.data.frame(x)
    aux = selection(x, y, q = qT, prevar = prevars, method = optionT,
                    family = family, seconds = FALSE, criterion = "deviance", nfolds = 2, cluster = FALSE)
    # if (!exists("pred")) {pred <<- aux$Prediction}
    pred <- aux$Prediction
    sel_num <- aux$Variable_number

    #res = y - pred
    res = aux$Best_model$residuals
    if (speed == TRUE & qT!=(nvar-1)) {
      xno = x[, -sel_num]
      var_imp = selection(xno, res, q = 1, method = optionT,
                          family = "gaussian", seconds = FALSE,
                          criterion = "deviance", cluster = FALSE)$Variable_number
      xres = xno[, c(var_imp)]
      }else if (speed == TRUE & qT == (nvar-1)){
       xres = x[, -sel_num]
      } else {
          xno = x[, -sel_num]
          realqmin <- qmin-qT
          var_imp = selection(xno, res, q = realqmin, method = optionT,
                             family = "gaussian", seconds = FALSE,
                            criterion = "deviance", cluster = FALSE)$Variable_number
          xres = xno[, c(var_imp)]
        }
    data_res = cbind(res, xres)



    if(ncol(data_res) == 2){
      if (optionT == "gam" & is.factor(xres) == FALSE) {
        xnam = paste("s(xres)", sep = "")
      } else {
        xnam = paste("xres", sep = "")
      }
    }else{
      xnam <- c()
      for (num in 1:(ncol(data_res)-1)) {
        if (optionT == "gam" & is.factor(data_res[, num+1]) == FALSE) {
          xnam[num] = paste("s(xres[,", num, "])", sep = "")
        } else {
          xnam[num] = paste("xres[,", num, "]", sep = "")
        }
      }
    }



    if(optionT == "gam"){
      fmla <- as.formula(paste("res ~ ", paste(xnam, collapse= "+")))
      pred1 <- gam(fmla, data = as.data.frame(data_res))
    }else{
      pred1 <- glm(res ~ ., family = "gaussian", data = as.data.frame(data_res))
    }

    pred1 <- predict(pred1, type = "response")
    T = sum(abs(pred1))
   # print(T)
    return(list(T = T, pred = pred, sel_num = sel_num))
  }




  if (missing(x)) {
    stop("Argument \"x\" is missing, with no default")
  }
  if (missing(y)) {
    stop("Argument \"y\" is missing, with no default")
  }
  if (speedup == "FALSE" & is.null(qmin)) {
    stop("Argument \"qmin\" is missing, with no default")
  }
  if (speedup == "TRUE" & !is.null(qmin)) {
    warning("Argument \"qmin\" has not been considered")
  }
  if (unique == "TRUE" & is.null(q)) {
    stop("Argument \"q\" is missing, with no default")
  }



  if (cluster == TRUE & detectCores() == 2 & is.null(ncores)) {
    stop("The number of cores used in the parallelized procedure is just one. It is recommended to use cluster = FALSE ")
  }

  # for paralellize
  if (cluster == TRUE){
    if (is.null(ncores)){
      ncores <- detectCores() - 1
    }else{
      ncores <- ncores
    }
    if(.Platform$OS.type == "unix"){par_type = "FORK"}else{par_type = "PSOCK"}
    cl <- makeCluster(ncores, type = par_type)
    on.exit(stopCluster(cl))
  }

  nvar <- ncol(x)
  n <- length(y)
  xydata <- cbind(y, x)
  pvalue <- c()
  Decision <- c()
  Hypothesis <- c()
  T <- c()
  ii <- 1
  if (unique == FALSE) {
    bucle <- c(1:(nvar - 1))
  } else {
    bucle <- q
  }
  pre <- c()
  for (qh0 in bucle) {
    print(paste("Processing IC bootstrap for H_0 (",
                qh0, ")..."), sep = "")
    if(isTRUE(unique)){pre == NULL}else{
      if(qh0 == 1){pre = NULL}else{pre = sel_numg}
    }
    #pred <- c()
    #sel_num <- c()
    res <- Tvalue(xy = xydata, qT = qh0, prevars = pre)
    T[ii] <- res$T
    sel_numg<- res$sel_num # lo saco de la funcion Tvalue bajo H_0
    muhatg <- res$pred  #lo saco de la funcion Tvalue bajo H_0



    if(family != "gaussian") muhatg[muhatg < 0] <- 0  #VER ESTO!!!!


    # Bootstrap
    ################################

    if (!is.null(bootseed)) {set.seed(bootseed)}

    if (family == "gaussian") {
      errg = y - muhatg
      err1 = errg * (1 - sqrt(5))/2
      err2 = errg * (1 + sqrt(5))/2

      funreplicate <- function(){
        yaux <- rbinom(n, 1, prob = (5 + sqrt(5))/10)
        return(muhatg + (err1 * yaux + err2 * (1 - yaux)))
      }

      yb <- replicate(nboot,funreplicate())
    }

    if (family == "binomial") {
      funreplicatebinom <- function(){
        yaux <- rbinom(n, 1, prob = muhatg)
        return(yaux)
      }
      yb <- replicate(nboot,funreplicatebinom())
    }

    if (family == "poisson") {
      funreplicatepois <- function(){
        yaux <- rpois(n, lambda = muhatg)
        return(yaux)
      }
      yb <- replicate(nboot,funreplicatepois())
    }

    funapply<-function(y){Tvalue(xy = cbind(y,xydata[,-1]), qT = qh0, prevars = NULL)$T}

    if (cluster == TRUE){
      Tboot <- parCapply(cl=cl,yb,funapply)
    }else{
      Tboot <- apply(yb,2,funapply)
    }

    pvalue[ii] = sum(Tboot > T[ii])/nboot

    if (pvalue[ii] >= 0.05) {
      Decision[ii] = "Not Rejected"
    } else {
      Decision[ii] = "Rejected"
    }
    Hypothesis[ii] = paste("H_0 (", qh0, ")", sep = "")
    T[ii] = round(T[ii], 2)
    ii = ii + 1
    if (Decision[ii - 1] == "Not Rejected") {break}
  }

  m = cbind(Hypothesis = Hypothesis, Statistic = T, pvalue = pvalue, Decision = Decision)
  cat("\n*************************************\n")
  res = as.data.frame(m)
  nvar = length(pvalue)
  print(res)
  res2 <- list(table = res, nvar = nvar)
  #return(res2)
}
