#' Add a nonlinear function to the model
#' 
#' This function takes a model and adds a non-linear function if
#' the likelihood-ratio supports this (via the 
#' \code{\link[stats]{anova}(..., test="chisq")} test for \pkg{stats}
#' while for \pkg{rms} you need to use the \code{rcs} spline
#' that is automatically evaluated for non-linearity).
#' 
#' @param model The model that is to be evaluated and adatpted for non-linearity
#' @param variable The name of the parameter that is to be tested for non-linearity.
#'  \emph{Note} that the variable should be included plain (i.e. as a linear variable)
#'  form in the model.
#' @param spline_fn Either a string or a function that is to be used
#'  for testing alternative non-linearity models
#' @param flex_param A \code{vector} with values that are to be tested as the 
#'  default second parameter for the non-linearity function that you want to 
#'  evaluate. This defaults to 2:7, for the \code{\link[splines]{ns}} it tests 
#'  the degrees of freedom ranging between 2 and 7.
#' @param min_fn This is the function that we want to minmized if the variable supports
#'  the non-linearity assumption. E.g. \code{\link[stats]{BIC}} or 
#'  \code{\link[stats]{AIC}}, note that the \code{BIC} will in the majority of cases
#'  support a lower complexity than the \code{AIC}.
#' @param sig_level The significance level for which the non-linearity is deemed
#'  as significant, defaults to 0.05.
#' @param verbal Set this to \code{TRUE} if you want print statements with the
#'  anova test and the chosen knots.
#' @param workers The function tries to run everything in parallel. Under some
#'  circumstances you may want to restrict the number of parallel threads to less
#'  than the defaul \code{dectectCores() - 1}, e.g. you may run out of memory
#'  then you can provide this parameter. If you do not want to use parallel then 
#'  simply set workers to \code{FALSE}. The cluster created using \code{\link[parallel]{makeCluster}} 
#'  function.
#' @param ... Passed onto internal \code{\link{prNlChooseDf}} function.
#' 
#' @example inst/examples/addNonlinearity_example.R
#' @rdname addNonlinearity
#' @import stats
#' @export
addNonlinearity <- 
  function(model, 
           variable,
           spline_fn,
           flex_param = 2:7,
           min_fn = AIC, 
           sig_level = .05,
           verbal = FALSE,
           workers,
           ...){
    UseMethod("addNonlinearity")  
  }

#' @export
addNonlinearity.default <- function(model, 
                                    variable,
                                    spline_fn,
                                    flex_param = 2:7,
                                    min_fn = AIC, 
                                    sig_level = .05,
                                    verbal = FALSE,
                                    workers,
                                    ...){
  stop("Model of class '", paste(class(model),
                                 collapse = ">"), "'",
       " not yeat implemented. Feel free to check and send me a suggestion.")
}

#' @export
addNonlinearity.lm <- function(model, 
                               ...){
  # Works fine with the glm method
  addNonlinearity.glm(model, ...)
}

#' @rdname addNonlinearity
#' @export
addNonlinearity.negbin <- function(model, 
                                   ...){
  # Works fine with the glm method
  addNonlinearity.glm(model, ..., libraries = "MASS")
}

#' @export
addNonlinearity.glm <- function(model, 
                                variable,
                                spline_fn,
                                flex_param = 2:7,
                                min_fn = AIC, 
                                sig_level = .05,
                                verbal = FALSE,
                                workers,
                                ...)
{
  if (is.function(spline_fn))
    spline_fn <- deparse(spline_fn)
  
  simplest_nonlinear <-
    eval(update(model,
                sprintf(".~.-%s+%s(%s, %d)",
                        variable, spline_fn, variable, min(flex_param)),
                x=FALSE, y=FALSE, evaluate = FALSE),
         envir = environment(formula(model)))
  anova_rslt <- anova(model,
                      simplest_nonlinear,
                      test = "Chisq")
  
  if  (verbal){
    cat("\n * Non-linearity for", variable, "*\n")
    print(anova_rslt)
  }
  
  # No evidence for non-linearity
  if (tail(anova_rslt[grep("Pr(", names(anova_rslt), fixed=TRUE)], 1) > sig_level){
    return(model)
  }
  
  return(prNlChooseDf(model = model, 
                      flex_param = flex_param, 
                      variable = variable, 
                      spline_fn = spline_fn, 
                      min_fn = min_fn, 
                      simplest_nonlinear = simplest_nonlinear, 
                      verbal = verbal,
                      workers = workers,
                      ...))
}

#' @export
addNonlinearity.coxph <- function(model, 
                                  variable,
                                  spline_fn,
                                  flex_param = 2:7,
                                  min_fn = AIC, 
                                  sig_level = .05,
                                  verbal = FALSE,
                                  workers,
                                  ...){
  if (is.function(spline_fn))
    spline_fn <- deparse(spline_fn)
  
  if (spline_fn == "pspline"){
    nonlinear <-
      eval(update(model,
                  sprintf(".~.-%s+%s(%s)",
                          variable, spline_fn, variable),
                  x=FALSE, y=FALSE, evaluate = FALSE),
           envir = environment(formula(model)))
    anova_rslt <- anova(model,
                        nonlinear,
                        test = "Chisq")
  }else{
    simplest_nonlinear <-
      eval(update(model,
                  sprintf(".~.-%s+%s(%s, %d)",
                          variable, spline_fn, variable, min(flex_param)),
                  x=FALSE, y=FALSE, evaluate = FALSE),
           envir = environment(formula(model)))
    anova_rslt <- anova(model,
                        simplest_nonlinear,
                        test = "Chisq")
  }
  if  (verbal){
    cat("\n * Non-linearity for", variable, "*\n")
    print(anova_rslt)
  }
  
  # No evidence for non-linearity
  if (tail(anova_rslt[grep("P(", names(anova_rslt), fixed=TRUE)], 1) > sig_level){
    if  (verbal){
      cat("\nNo support for nonlinearity for variable: '", variable, "' as defined by the criteria",
          sep = "")
    }
    return(model)
  }
  
  if (spline_fn == "pspline")
    return(nonlinear)
  
  return(prNlChooseDf(model = model, 
                      flex_param = flex_param, 
                      variable = variable, 
                      spline_fn = spline_fn, 
                      min_fn = min_fn, 
                      simplest_nonlinear = simplest_nonlinear, 
                      verbal = verbal,
                      workers = workers,
                      libraries = "survival"))
}

#' @export
addNonlinearity.rms <-
  function(model, 
           variable,
           spline_fn,
           flex_param = 2:7,
           min_fn = AIC, 
           sig_level = .05,
           verbal = FALSE,
           workers,
           ...){
    if (is.function(spline_fn))
      spline_fn <- deparse(spline_fn)
    
    if (spline_fn != "rcs")
      stop("Only works with the rcs() function of the rms-package")
    
    simplest_nonlinear <-
      eval(update(model,
                  sprintf(".~.-%s+%s(%s, %d)",
                          variable, spline_fn, variable, min(flex_param)),
                  x=FALSE, y=FALSE, evaluate = FALSE),
           envir = environment(formula(model)))
    
    anova_rslt <- anova(simplest_nonlinear)
    
    row <- which(variable == substr(rownames(anova_rslt), 1, nchar(variable)))
    if (!grepl("Nonlinear", rownames(anova_rslt)[row + 1]))
      stop("Expect a 'Nonlinear' row after the spline row")
    
    if  (verbal){
      cat("\n * Non-linearity for", variable, "*\n")
      print(anova_rslt[0:1 + row, ])
    }
    
    if (anova_rslt[row+1,"P"] > sig_level){
      return(model)
    }
    
    return(prNlChooseDf(model = model, 
                        flex_param = flex_param, 
                        variable = variable, 
                        spline_fn = spline_fn, 
                        min_fn = min_fn, 
                        simplest_nonlinear = simplest_nonlinear, 
                        verbal = verbal,
                        workers = workers,
                        libraries = "rms"))
  }

#' Chooses the degrees of freedom for the non-linearity
#' 
#' Looks for the model with the minimal \code{min_fn} within the flex_param span.
#' 
#' @param libraries If we use the parallel approach we need to make sure that the 
#'  right libraries are available in the threads
#' @param simplest_nonlinear The simplest non-linear form that the ANOVA has been tested against
#' @import parallel
#' @inheritParams addNonlinearity
#' @keywords internal
prNlChooseDf <- function (model, 
                          flex_param, 
                          variable, 
                          spline_fn, 
                          min_fn, 
                          simplest_nonlinear, 
                          verbal,
                          workers,
                          libraries) {
  local_workers <- FALSE
  if (missing(workers)){
    local_workers <- TRUE
    
    # Look for best BIC using a parallel approach
    workers <- makeCluster(max(detectCores()-1, 1))
    
    if (!is.null(model$call$data))
      clusterExport(workers, deparse(model$call$data),
                    envir = environment(formula(model)))
    
    
    tmp <- clusterEvalQ(workers, library(splines))
    
    # Load libraries necessary into the workers
    if (!missing(libraries)){
      for (l in libraries){
        tmp <- sprintf("clusterEvalQ(workers, library(%s))", l) %>% 
          parse(text = .) %>% 
          eval()
      }
    }
  }
  
  if (is.logical(workers) && 
      workers == FALSE){
    fits <-
      lapply(X = flex_param[-which.min(flex_param)],
             FUN = 
               function(x, model, variable, spline_fn) {
                 eval(update(model,
                             sprintf(".~.-%s+%s(%s, %d)",
                                     variable, spline_fn, variable, x),
                             y = FALSE, x = FALSE, evaluate = FALSE),
                      envir = environment(formula(model)))
               },
             model = model,
             variable = variable,
             spline_fn = spline_fn)
  }else{
    fits <-
      parLapply(cl = workers,
                X = flex_param[-which.min(flex_param)],
                function(x, model, variable, spline_fn) {
                  eval(update(model,
                              sprintf(".~.-%s+%s(%s, %d)",
                                      variable, spline_fn, variable, x),
                              y = FALSE, x = FALSE, evaluate = FALSE),
                       envir = environment(formula(model)))
                },
                model = model,
                variable = variable,
                spline_fn = spline_fn)
    if (local_workers){
      stopCluster(workers)
    }
  }
  
  
  if (!is.function(min_fn))
    min_fn <- get(min_fn)
  
  # Add the simplest form
  fits <- c(list(simplest_nonlinear), fits)
  
  no_knots <- which.min(sapply(fits, min_fn))
  if  (verbal){
    if (no_knots == 1){
      cat("\n -- Chose the simplest form with", min(flex_param), "flex\n")
    }else{
      cat("\n -- Chose:", flex_param[-which.min(flex_param)][no_knots-1], "flex\n")
    }
  }
  
  return(fits[[no_knots]])
}