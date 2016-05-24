##' Eval a single list of specifications
##' 
##' Evaluate a single list of specifications. Climate parameters are selected
##' from the parameter matrix by name, then the appropriate function is applied
##' to the correct month subset.
##' @param clim climate data as returned from make_pmat
##' @param sel single list of parameter specifications
##' @return a list
##' @keywords internal manip
eval_selection <- function(clim, sel) {

  aggregate <- list()
  NAMES <- character(0)

  ## methods for parameter aggregation using .method, .month and .param
  method_full <- function(x) {
    x
  }

  method_mean <- function(x) {
    rowMeans(x)
  }

  method_sum <- function(x) {
    rowSums(x)
  }
  
  n <- length(as.list(sel))

  is_method_argument <- function(x) {
    any(x[1] == c("mean", "sum", "full"))
  }
  
  is_param_argument <- function(x) {
    if (is.character(x[1]) && !any(x[1] == c("mean", "sum", "full"))) {
      if (any(x[1] == names(clim))) {
        TRUE
      } else {
        stop(paste("Parameter name '", x[1], "' not recognized.", sep
                   = ""))
      }
    } else {
      FALSE
    } 
  }

  ## default values
  month <- format_month(-6:9)
  param <- NULL
  method <- "full"
  if (n > 0) {
    if (any(sapply(sel, is.numeric))) {
      ## get numeric argument
      month <- format_month(unlist(sel[sapply(sel, is.numeric)]))
    }
    if (any(sapply(sel, is_method_argument))) {
      ## method as argument
      method <- unlist(sel[sapply(sel, is_method_argument)])
    }
    if (any(sapply(sel, is_param_argument))) {
      ## param as argument
      param <- unlist(sel[sapply(sel, is_param_argument)])
    }
  }

  ## param == NULL means that all parameters are used
  if (is.null(param)) {
    param <- names(clim)
  }

  for (i in seq_len(length(param))) {   # do this for all climate
    ## parameters

    ## select climate element by climate parameter name
    eval(
      substitute(
        this_clim <- clim$name,
        list(name = param[i])
        )
      )
    
    ## and then select only the relevant months
    selected_clim <- this_clim[,month$match]

    ## switch on method and apply to climate
    .method <- switch(method,
                      full = method_full,
                      mean = method_mean,
                      sum = method_sum)
    
    .clim <- .method(selected_clim)

    ## generate variable names
    .names <- switch(method,
                     full = paste(param[i], month$names, sep = "."),
                     mean = paste(param[i], "mean", paste(month$names,
                       collapse = "."), sep = "."),
                     sum = paste(param[i], "sum", paste(month$names,
                       collapse = "."), sep = "."),)
    aggregate[[i]] <- .clim
    NAMES <- c(NAMES, .names)
  }

  ## combine all climate parameters into a flat data.frame
  aggregate <- as.data.frame(aggregate)

  ## prepare return value
  ret <- list()
  ret$month <- month                    # the months used
  ret$method <- method                  # the method used
  ret$param <- param                    # the parameter used
  ret$aggregate <- aggregate            # the data itself
  ret$names <- NAMES                    # the names of the aggregated
  ## variables
  ret
}
