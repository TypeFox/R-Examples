#' Compares different scores in different regression objects.
#' 
#' Creates a composite from different regression objects into
#' one forestplot where you can choose the variables of interest
#' to get an overview and easier comparison.
#' 
#' @param regr.obj A list with all the fits that have variables that are to 
#'   be identified through the regular expression 
#' @param variablesOfInterest.regexp A regular expression identifying the variables
#'   that are of interest of comparing. For instance it can be "(score|index|measure)"
#'   that finds scores in different models that should be compared.
#' @param reference.names Additional reference names to be added to each model
#' @param rowname.fn A function that takes a rowname and sees if it needs
#'   beautifying. The function has only one parameter the coefficients name and should
#'   return a string or expression.  
#' @param estimate.txt The text of the estimate, usually HR for hazard ratio, OR for 
#'  odds ratio
#' @param zero Indicates what is zero effect. For survival/logistic fits the zero is
#'   1 while in most other cases it's 0.
#' @param exp Report in exponential form. Default true since the function was built for 
#'   use with survival models.
#' @param add_first_as_ref If you want that the first variable should be reference for 
#'   that group of variables. The ref is a variable with the estimate 1 or 0 depending
#'   if exp() and the confidence interval 0.
#' @param ref_labels If add_first_as_ref is TRUE then this vector is used for the model
#'   fits. 
#' @param ref_txt Text instead of estimate number
#' @param digits Number of digits to use for the estimate output
#' @param ... Passed to \code{\link[forestplot]{forestplot}}
#' 
#' @example inst/examples/forestplotCombineRegrObj_example.R
#' 
#' @inheritParams forestplot::forestplot
#' @import forestplot
#' @family \code{\link[forestplot]{forestplot}} wrappers
#' @export
forestplotCombineRegrObj <- function( 
  regr.obj,
  variablesOfInterest.regexp,
  reference.names,
  rowname.fn,
  estimate.txt,
  exp              = xlog,
  add_first_as_ref = FALSE,
  ref_txt          = "ref.",
  ref_labels       = c(),
  digits           = 1,
  is.summary,
  xlab,
  zero,
  xlog,
  ...)
{   
  if (length(regr.obj) < 2)
    stop("This function combines several fits so please provide more than one model_fits")
  
  if (is.object(regr.obj) == TRUE)
    stop("The model_fits need to be a list of fits")

  if (!missing(reference.names) &&
        length(reference.names) != length(regr.obj))
    stop("The number or reference names need to be equal to the number of",
         " regression model fits. You have ", length(reference.names), " references ",
         " and ", length(regr.obj), " regression objects")
    
  # Initiate some standard values if the user
  # hasn't supplied any
  if (missing(xlab)){
    if (isFitCoxPH(regr.obj[[1]])){
      xlab <- "Hazard Ratio"
    }else if(isFitLogit(regr.obj[[1]])){
      xlab <- "Odds Ratio"
    }
  }
  
  if (missing(estimate.txt)){
    if (isFitCoxPH(regr.obj[[1]])){
      estimate.txt <- "HR"
    }else if(isFitLogit(regr.obj[[1]])){
      estimate.txt <- "OR"
    }
  }
  
  if (missing(xlog)){
    if (isFitCoxPH(regr.obj[[1]]) ||
          isFitLogit(regr.obj[[1]])){
      xlog <- TRUE
      if (missing(zero))
        zero <- 1
      if (missing(exp))
        exp <- TRUE
    }else{
      xlog <- FALSE
      if (missing(zero))
        zero <- 0
      if (missing(exp))
        exp <- FALSE
    }
  }
  
  models_fit_fp_data = getModelData4Forestplot(regr.obj = regr.obj,
                                               exp = exp,
                                               variablesOfInterest.regexp = variablesOfInterest.regexp,
                                               ref_labels = ref_labels,
                                               add_first_as_ref = add_first_as_ref)
  
  rn <- list()
  t.coef <- c(NA)
  t.low <- c(NA)
  t.high <- c(NA)
  t.is_ref <- c(NA)
  for(i in 1:length(models_fit_fp_data)){
    if (length(rn) > 0){
      rn <- append(rn, NA)
      t.coef <- append(t.coef, NA)
      t.low <- append(t.low, NA)
      t.high <- append(t.high, NA)
      t.is_ref <- append(t.is_ref, FALSE)
    }
    if (!missing(reference.names)){
      rn <- append(rn, reference.names[i])
      t.coef <- append(t.coef, NA)
      t.low <- append(t.low, NA)
      t.high <- append(t.high, NA)
    }
    raw_row_names <- as.list(rownames(models_fit_fp_data[[i]]))
    if (!missing(rowname.fn)){
      if (is.character(rowname.fn))
        rowname.fn <- get(rowname.fn)
      
      for(c in 1:length(raw_row_names))
        raw_row_names[[c]] <- rowname.fn(raw_row_names[[c]])
    }
    rn <- append(rn, raw_row_names)
    t.coef <- append(t.coef, models_fit_fp_data[[i]][, "beta"])
    t.low <- append(t.low, models_fit_fp_data[[i]][, "low"])
    t.high <- append(t.high, models_fit_fp_data[[i]][, "high"])
    t.is_ref <- append(t.is_ref, is.na(models_fit_fp_data[[i]][, "p_val"]))
  }
  
  rn <- append("Variable", rn)
  
  col2 <- list(estimate.txt)
  # TODO: should probably use options("digits") but
  # it defaults to 7 why this might be more annoying
  # than helpful
  for (i in 2:length(t.coef))
    col2 <- append(col2, ifelse(is.na(t.coef[i]), "", sprintf(sprintf(" %%.%df ", digits), t.coef[i])))
    
  if (add_first_as_ref)
    col2[t.is_ref] <- ref_txt
  
  rn <- list(
    rn,
    col2)
  
  if (missing(is.summary))
    is.summary <- c(TRUE, rep(FALSE, length(t.coef)-1))
  
  forestplot(rn, 
    xlim       = c(0,10),
    mean       = t.coef, 
    lower      = t.low, 
    upper      = t.high,
    is.summary = is.summary,
    zero       = zero,
    xlab       = xlab,
    xlog       = xlog,
    ...)
  
}
