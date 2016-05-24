lm.nri <- function(formula, preddata=NULL, ...)
{
  lm_apply_lh <- function(response, ...)
  {
    formula <- get("formula", envir= environment_apply)
    predictors <- get("predictors", envir= environment_apply)
    data <- cbind(data.frame(response=response), predictors)
    model <- lm(formula = formula, data = data, ...)
    r.squared <- get("r.squared", envir= environment_apply)
    r.squared <- c(r.squared, summary(model)$r.squared)
    assign("r.squared", r.squared, environment_apply)
    return(summary(model)$coefficients[,c(1:4)])
  }
  lm_apply_rh <- function(predictor, ...)
  {
    formula <- get("formula", envir= environment_apply)
    response <- get("response", envir= environment_apply)
    data <- data.frame(response=response[,1], predictor=predictor)
    model <- lm(response ~ predictor, data = data, ...)
    su <- summary(model)$coefficients
    if (nrow(su)<ncol(data))
      su <- rbind(su, rep.int(NA,4))
    r.squared <- get("r.squared", envir= environment_apply)
    r.squared <- c(r.squared, summary(model)$r.squared)
    assign("r.squared", r.squared, environment_apply)
    return(su[,c(1:4)])
  }
  nri_response <- 0
  mc <- match.call()

  formula <- as.formula(formula)
  mc <- formula



  if (mode(preddata)!="NULL")
  {
    if (is.speclib(preddata))
      preddata <- preddata$attributes
    data <- preddata
    formula <- terms(formula, data=data)
  } else {
    formula <- terms(formula)
    data <- NULL
  }
  vars <- all.vars(formula)
  for (i in 1:length(vars))
  {
    x <- try(eval(parse(text = vars[i])),
             silent = TRUE
            )
    if (class(x) == "try-error")
      x <- eval(parse(text = vars[i]), data, environment(formula))
    if (class(x)=="Nri")
    {
      if (is.element(vars[i], attr(formula,"term.labels")))
      {
        nri_response <- 1
      } else {
        if (!is.element(vars[i],row.names(attr(formula,"factors"))))
          stop("Could not determine which variable contains nri-values. This may be caused by a function to be applied to nri-values")
        nri_response <- -1
      }
    }
  }

  if (nri_response==0)
    stop("Could not determine which variable contains nri-values")

  if (nri_response==-1)
  {
    formula_apply <- update.formula(mc, response ~ .)
    response <- as.character(mc)[2]
    x <- eval(parse(text = response))
    response <- paste(response,"nri",sep="$")
    response <- eval(parse(text = response))

    mat <- match(all.vars(formula_apply)[-1], names(data))

    if (length(mat)==1)
    {
      predictors <- data.frame(a=data[,mat])
      names(predictors) <- names(data)[mat]
    } else {
      predictors <- data[,mat]
    }

    environment_apply <- new.env(parent = .GlobalEnv)

    assign("formula", formula_apply, environment_apply)
    assign("predictors", predictors, environment_apply)
    assign("r.squared", vector(mode="numeric", length=0), environment_apply)

    lm_data <- apply(response, MARGIN = c(1, 2), FUN = lm_apply_lh, ...)

  } else {
    if (length(vars) > 2)
      stop("Models with more than 2 variables with nri values as predictors not supported, yet")
    predictor <- as.character(mc)[3]
    x <- eval(parse(text = predictor))
    predictor <- paste(predictor,"nri",sep="$")
    predictor <- eval(parse(text = predictor))

    formula_apply <- update.formula(mc, . ~ predictor)

    mat <- match(all.vars(formula_apply)[1], names(data))

    response <- data.frame(a=data[,mat])
    names(response) <- names(data)[mat]

    environment_apply <- new.env(parent = .GlobalEnv)

    assign("formula", formula_apply, environment_apply)
    assign("response", response, environment_apply)
    assign("r.squared", vector(mode="numeric", length=0), environment_apply)

    lm_data <- apply(predictor, MARGIN = c(1, 2), FUN = lm_apply_rh, ...)
  }

  ncol = dim(lm_data)[1]/4

  nam <- NULL
  nam[[1]] <- c("(Intercept)", attr(formula,"term.labels"))
  nam[[2]] <- paste("B_",x$wavelength, sep="")
  nam[[3]] <- paste("B_",x$wavelength, sep="")
  final <- list(estimate = distMat3D(as.numeric(t(lm_data[1:ncol,])),
                                     length(x$wavelength),
                                     ncol),
                std.error = distMat3D(as.numeric(t(lm_data[c((ncol+1):(2*ncol)),])),
                                     length(x$wavelength),
                                     ncol),
                t.value = distMat3D(as.numeric(t(lm_data[c((2*ncol+1):(3*ncol)),])),
                                     length(x$wavelength),
                                     ncol),
                p.value = distMat3D(as.numeric(t(lm_data[c((3*ncol+1):(4*ncol)),])),
                                     length(x$wavelength),
                                     ncol),
                r.squared = distMat3D(as.numeric(get("r.squared",
                                                     envir= environment_apply)),
                                      ncol = length(x$wavelength),
                                      nlyr = 1)
           )
  attr(final,"call") <- mc
  attr(final,"function") <- "lm"
  if (nri_response==-1)
  {
    attr(final,"is.predictor.nri") <- FALSE
    attr(final,"predictors") <- predictors
  } else {
    attr(final,"is.predictor.nri") <- TRUE
    attr(final,"response") <- response
  }
  x@multivariate <- final
  return(x)
}
