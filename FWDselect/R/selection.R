#'Selecting a subset of \code{q} variables
#'
#'@description Main function for selecting the best subset of \eqn{q} variables.
#'  Note that the selection procedure can be used with lm, glm or gam functions.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param q An integer specifying the size of the subset of variables to be
#'  selected.
#'@param prevar A vector containing the number of the best subset of
#'  \code{q-1} variables. \code{NULL}, by default.
#'@param criterion The information criterion to be used.
#'  Default is the deviance. Other functions provided
#'  are the coefficient of determination (\code{"R2"}), the residual
#'  variance (\code{"variance"}), the Akaike information criterion (\code{"aic"}),
#'  AIC with a correction for finite sample sizes (\code{"aicc"})
#'  and the Bayesian information criterion (\code{"bic"}). The deviance,
#'  coefficient of determination and variance are calculated by cross-validation.
#'@param method A character string specifying which regression method is used,
#'  i.e., linear models (\code{"lm"}), generalized additive models
#'  (\code{"glm"}) or generalized additive models (\code{"gam"}).
#'@param family A description of the error distribution and link function to be
#'  used in the model: (\code{"gaussian"}), (\code{"binomial"}) or
#'  (\code{"poisson"}).
#'@param seconds A logical value. By default, \code{FALSE}.  If \code{TRUE}
#'  then, rather than returning the single best model only, the function returns
#'  a few of the best models (equivalent).
#'@param nmodels Number of secondary models to be returned.
#'@param nfolds Number of folds for the cross-validation procedure, for
#'\code{deviance}, \code{R2} or \code{variance} criterion.
#'@param cluster A logical value. If  \code{TRUE} (default), the
#'  procedure is  parallelized. Note that there are cases without enough
#'  repetitions (e.g., a low number of initial variables) that R will gain in
#'  performance through serial computation. R takes time to distribute tasks
#'  across the processors also it will need time for binding them all together
#'  later on. Therefore, if the time for distributing and gathering pieces
#'  together is greater than the time need for single-thread computing, it does
#'  not worth parallelize.
#'@param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores to be used
#' is equal to the number of cores of the machine - 1.
#'@return
#'\item{Best model}{The best model. If \code{seconds=TRUE}, it returns
#'  also the best alternative models.}
#'  \item{Variable name}{Names of the variable.}
#'  \item{Variable number}{Number of the variables.}
#'  \item{Information criterion}{Information criterion used and its value.}
#'  \item{Prediction}{The prediction of the best model.}
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' obj1 = selection(x, y, q = 1, method = "lm", criterion = "variance", cluster = FALSE)
#' obj1
#'
#' # second models
#' obj11 = selection(x, y, q = 1, method = "lm", criterion = "variance",
#' seconds = TRUE, nmodels = 2, cluster = FALSE)
#' obj11
#'
#' # prevar argument
#' obj2 = selection(x, y, q = 2, method = "lm", criterion = "variance", cluster = FALSE)
#' obj2
#' obj3 = selection(x, y, q = 3, prevar = obj2$Variable_numbers,
#' method = "lm", criterion = "variance", cluster = FALSE)
#'
#'
#'@importFrom mgcv gam
#'@importFrom mgcv predict.gam
#'@importFrom parallel detectCores
#'@importFrom parallel makeCluster
#'@importFrom parallel parLapply
#'@importFrom parallel stopCluster
#'@importFrom stats as.formula
#'@importFrom stats deviance
#'@importFrom stats lm
#'@importFrom stats glm
#'@importFrom stats predict
#'@importFrom stats update
#'@importFrom stats var
#'@importFrom stats AIC
#'@importFrom stats BIC
#'@importFrom stats logLik
#'@export


selection <- function(x, y, q, prevar = NULL, criterion = "deviance",
                      method = "lm", family = "gaussian", seconds = FALSE,
                      nmodels = 1, nfolds = 5, cluster = TRUE, ncores = NULL) {



  if (missing(x)) {
    stop("Argument \"x\" is missing, with no default")
  }
  if (missing(y)) {
    stop("Argument \"y\" is missing, with no default")
  }
  if (missing(q)) {
    stop("Argument \"q\" is missing, with no default")
  }

  nvar <- ncol(x)
  inside <- integer(q)
  n <- length(y)


  if(q == nvar) {
    stop('The size of subset \'q\' is the same that the number of covariates')
  }


  if(!criterion %in% c("deviance", "R2", "variance", "aic", "aicc", "bic")) {
    stop('The selected criterion is not implemented')
  }


  if (cluster == TRUE & detectCores() == 2 & is.null(ncores)) {
    stop("The number of cores used in the parallelized procedure is just one.
         It is recommended to use cluster = FALSE ")
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

  #dat = data.frame(y,x)

  if (method == "lm") {
    model <- lm(y ~ NULL)
  }
  if (method == "glm") {
    model <- glm(y ~ NULL, family = family)
  }
  if (method == "gam") {
    model <- gam(y ~ NULL, family = family)
  }



  # To use the variable of the previous q and
  # it have not to look for again (class(prevar) = vector)
  if (is.null(prevar)) { }else{
    xyes = c()
    for (l in 1:(q-1)){
      if (method == "gam" & is.factor(x[, prevar[l]]) == FALSE) {
        xnam = paste("s(x[,", prevar[l], "])", sep = "")
      } else {
        xnam = paste("x[,", prevar[l], "]", sep = "")
      }
      xyes[l] = xnam
    }

    form1 <- update(as.formula(model, env = environment(fun = NULL)), paste(". ~ ",
                                                                            paste(xyes, collapse = "+")))
    if (method == "gam"){
      model <- gam(form1, family = family)
    }else{
      model <- glm(form1, family = family)
    }

    # model <- update(model, as.formula(paste(". ~ ", paste(xyes, collapse = "+"))))

  }




  fwdstep <- function(j){

    form0 <- as.formula(model, env = environment(fun = NULL))

    if (method == "gam" & is.factor(x[,j]) == FALSE) {
      form1 <- update(form0, . ~ . + s(x[,j]))
    }else{
      form1 <- update(form0, . ~ . + x[,j])
    }


    if (method == "gam"){
      models <- gam(form1, family = family)
    }else{
      models <- glm(form1, family = family)
    }

    return(deviance(models))

  }




  fwdstep2 <- function(j, bucle){

    if (method == "gam" & is.factor(x[ ,j]) == FALSE) {
      xnam[bucle] <- paste("s(x[ ,", j, "])",sep="")
    } else {
      xnam[bucle] <- paste("x[ ,", j, "]",sep="")
    }

    form0 <- as.formula(model, env = environment())

    form1 <- update(form0, paste(". ~ ", paste(xnam, collapse = "+")))

    if (method == "gam"){
      model1 <- gam(form1, family = family)
    }else{
      model1 <- glm(form1, family = family)
    }

    return(deviance(model1))
  }



  out <- 1:nvar
  if(is.null(prevar)){
    xyes = NULL
    bucle <- c(1:q)
  }else{
    bucle <- q
    inside <- prevar
    out <- out[-prevar]
  }



  for (k in bucle) {
    ic <- NULL
    if (cluster == TRUE){
      ic <- parLapply(cl = cl, out, fwdstep)
    }else{
      ic <- sapply(out, fwdstep)
    }
    ii = which.min(ic)
    inside[k] = out[ii]
    out = out[-ii]

    if (method == "gam" & is.factor(x[, inside[[k]]]) == FALSE) {
      xnam = paste("s(x[,", inside[[k]], "])", sep = "")
    } else {
      xnam = paste("x[,", inside[[k]], "]", sep = "")
    }

    xyes[k] = xnam
    form1 <- update(as.formula(model, env = environment(fun = NULL)), paste(". ~ ",
                                                                            paste(xyes, collapse = "+")))
    if (method == "gam"){
      model <- gam(form1, family = family)
    }else{
      model <- glm(form1, family = family)
    }
    bestic = deviance(model)
  }


  ## Here it have introduced the first q variables



  stop <- integer(q)
  end <- 1
  if (q == 1 | q == nvar) {
    end <- 0
  }
  cont <- 0
  while (end != 0) {
    stop <- 0
    for (f in 1:q) {

      #para coger en un vector los nombres
      for (num in 1:length(inside)) {
        if (method == "gam" & is.factor(x[, inside[num]]) == FALSE) {
          xnam[num] = paste("s(x[,", inside[num], "])", sep = "")
        } else {
          xnam[num] = paste("x[,", inside[num], "]", sep = "")
        }
      }



      ic <- NULL
      if (cluster == TRUE){
        ic <- parLapply(cl = cl, out, fwdstep2, bucle = f)
      }else{
        ic <- sapply(out, fwdstep2, bucle = f)
      }

      ii = which.min(ic)
      if (ic[ii] >= bestic) {
        stop[f] <- 0
      } else {
        ii = which.min(ic)
        oldinside = inside
        inside[f] = out[ii]
        out[ii] = oldinside[f]
        if (method == "gam" & is.factor(x[ ,inside[f]]) == FALSE) {
          xin = paste("s(x[,", inside[f], "])", sep = "")
        } else {
          xin = paste("x[,", inside[f], "]", sep = "")
        }
        xnam[f] = xin
        #model <- update(model, as.formula(paste(". ~ ", paste(xnam, collapse = "+"))))
        form1 <- update(as.formula(model, env = environment()), paste(". ~ ", paste(xnam, collapse = "+")))
        if (method == "gam"){
          model <- gam(form1, family = family)
        }else{
          model <- glm(form1, family = family)
        }

        bestic = deviance(model)
        stop[f] = 1
      }
    }
    cont = cont + 1
    end = sum(stop)
  }

  pred <- predict(model, type = "response")


  # functions for cv
  cv <- function(nfolds){
    #function for calculate ic for each fold
    eachfold <- function(fold){
      test <- aux$which==fold
      Wtrainning = rep(1, n)
      Wtrainning[test] = 0
      formula <- eval(model$call$formula)
      dat <- data.frame(Wtrainning = Wtrainning)

      if (method == "lm") {
        Mtrainning = lm(formula, weights = Wtrainning, data = dat)
      }

      if (method == "glm") {
        Mtrainning = glm(formula, family = family, weights = Wtrainning, data = dat)
      }

      if (method == "gam") {
        Mtrainning = gam(formula, family = family, weights = Wtrainning, data = dat)
      }

      muhat = predict(Mtrainning, type = "response")
      muhat_test = muhat[test]
      y_test = y[test]
      if (family == "binomial") {y = as.numeric(as.character(y))}



      if (criterion == "deviance") {
        if (family == "gaussian"){
          dev_cv = sum((y_test - muhat_test)^2, na.rm = TRUE)
        }
        if (family == "binomial") {
          ii = muhat_test < 1e-04
          muhat_test[ii] = 1e-04
          ii = muhat_test > 0.9999
          muhat_test[ii] = 0.9999
          entrop = rep(0, length(test))
          ii = (1 - y_test) * y_test > 0
          if (sum(ii, na.rm = TRUE) > 0) {
            entrop[ii] = 2 * (y_test[ii] * log(y_test[ii])) +
              ((1 - y_test[ii]) * log(1 - y_test[ii]))
          } else {
            entrop = 0
          }
          entadd = 2 * y_test * log(muhat_test) +
            (1 - y_test) * log(1 - muhat_test)
          dev_cv = sum(entrop - entadd, na.rm = TRUE)
        }
        if (family == "poisson") {
          tempf = muhat_test
          ii = tempf < 1e-04
          tempf[ii] = 1e-04
          dev_cv = 2 * (-y_test * log(tempf) - (y_test - muhat_test))
          ii = y_test > 0
          dev_cv[ii] = dev_cv[ii] + (2 * y_test[ii] * log(y_test[ii]))
          dev_cv = sum(dev_cv, na.rm = TRUE)
        }
      } else if (criterion == "R2") {
        var_res = sum((y[test] - muhat[test])^2, na.rm = TRUE)/length(test)
        r2cv = 1 - (var_res/(var(y[test]) * (length(test) - 1)/length(test)))
      }else{
        var_res = sum((y[test] - muhat[test])^2, na.rm = TRUE)/length(test)
      }
      if (criterion == "deviance") {
        return(dev_cv)
      } else if (criterion == "R2") {
        return(r2cv)
      }else{
        return(var_res)
      }
    }

    aux <- cvTools::cvFolds(n, K = nfolds, type = "consecutive")
    if (cluster == TRUE){
      cv_ics <- parLapply(cl = cl, 1:nfolds, eachfold)
    }else{
      cv_ics <- sapply(1:nfolds, eachfold)
    }
    return(mean(unlist(cv_ics)))
  }

  aicc <- function(model){
    n <- length(model$y)
    k <- attr(logLik(model), "df")
    res <- AIC(model) + 2 * k * (k+1)/(n-k-1)
  }

  if(criterion %in% c("deviance", "R2", "variance")){
    icfin <- cv(nfolds)
  }else{
    if (criterion == "aic"){
     icfin <- AIC(model)
    }else if(criterion == "aicc"){
      icfin <- aicc(model)
    }else{
      icfin <- BIC(model)
    }
  }

  if(class(x) == "data.frame"){
    names1 = names(x[inside])
  }else{
    allnames <- colnames(x)
    names1 = allnames[inside]
  }

  if(is.null(names1)){names1=inside}  #por si no tiene nombres

  res <- list(Best_model = model, Variable_names = names1,
              Variable_numbers = inside, Information_Criterion = icfin,
              ic = criterion, seconds = seconds, nmodels = nmodels,
              Prediction = pred)



  # Second models

  if (seconds == TRUE) {
    bestic1 = bestic
    besticn = 0
    cont = -1
    fin = 1
    for (h in 1:nmodels) {
      cont = -1
      fin = 1
      while (fin != 0) {
        fin = 0
        for (zz in 1:q) {

          #para coger en un vector los nombres
          for (num in 1:length(inside)) {
            if (method == "gam" & is.factor(x[, inside[num]]) == FALSE) {
              xnam[num] = paste("s(x[,", inside[num], "])", sep = "")
            } else {
              xnam[num] = paste("x[,", inside[num], "]", sep = "")
            }
          }



          ic2 <- NULL
          if (cluster == TRUE){
            ic2 <- parLapply(cl = cl, out, fwdstep2, bucle = zz)
          }else{
            ic2 <- sapply(out, fwdstep2, bucle = zz)
          }
          ic2 <- unlist(ic2)


          if ((zz == 1) & (cont == -1)) {
            bestic = 1e+11}
          # oldinside = inside
          # inside[zz] = out[1]
          # out[1] = oldinside[1]
          # }

          for (j in 1:length(out)) {
            #  if ((zz == 1) & (cont == -1) &
            #       (j == 1)) {
            #   j = 2
            #  }
            if (h == 1) {
              if ((ic2[j] < bestic) &
                    (round(ic2[j],3) > round(bestic1,3))) {
                bestic = ic2[j]
                oldinside = inside
                inside[zz] = out[j]
                out[j] = oldinside[zz]
                fin = 1
              }
            } else {
              if ((ic2[j] < bestic) &
                    (ic2[j] > besticn)) {
                bestic = ic2[j]
                oldinside = inside
                inside[zz] = out[j]
                out[j] = oldinside[zz]
                fin = 1
              }
            }
          }
        }
        cont = cont + 1
      }


      for (num in 1:length(inside)) {
        if (method == "gam" & is.factor(x[, inside[num]]) == FALSE) {
          xnam[num] = paste("s(x[,", inside[num], "])", sep = "")
        } else {
          xnam[num] = paste("x[,", inside[num], "]", sep = "")
        }
      }

     # model <- update(model, as.formula(paste(". ~ ",paste(xnam, collapse = "+"))))

     form1 <- update(as.formula(model, env = environment()), paste(". ~ ", paste(xnam, collapse = "+")))
     if (method == "gam"){
       model <- gam(form1, family = family)
     }else{
       model <- glm(form1, family = family)
     }

      besticn = deviance(model)


      if(criterion %in% c("deviance", "R2", "variance")){
        icfin <- cv(nfolds)
      }else{
        if (criterion == "aic"){
          icfin <- AIC(model)
        }else if(criterion == "aicc"){
          icfin <- aicc(model)
        }else{
          icfin <- BIC(model)
        }
      }

      if(class(x) == "data.frame"){
        names2 = names(x[inside])
      }else{
        allnames <- colnames(x)
        names2 = allnames[inside]}
      if(is.null(names2)){names2=inside}  #por si no tiene nombres


      res2 <- list(Alternative_model = model, Variable_names = names2,
                   Variable_numbers = inside, Information_Criterion = icfin,
                   ic = criterion)


      res = c(res, res2)

    }
  }


  class(res) <- "selection"
  return(res)
}
