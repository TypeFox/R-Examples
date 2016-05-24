#' Get the confidence intervals
#' 
#' These are functions that get the estimates and the confidence intervals.
#' Due to package differences there are some local modifications.
#' 
#' @section The default:
#' 
#' Gets the estimate and confidence interval using the \code{\link[stats]{confint}}
#' and  \code{\link[stats]{coef}}.
#' 
#' @inheritParams getCrudeAndAdjustedModelData
#' @param skip_intercept If the model should remove the intercept from
#'  the returned values.
#'  
#' @return \code{matrix} Returns a n x 3 matrix where the n equals the number
#'  of variables.
#'  
#' @importFrom nlme intervals
#' @keywords internal
#' @rdname GetCoefAndCI
prCaDefaultGetCoefAndCI <- function(model, level, skip_intercept=FALSE){
  # Get the coefficients
  if (inherits(model, "lme")){
    tmp <- intervals(model, level=level)$fixed
    my_coefficients <- tmp[,"est."]
    coef_names <- rownames(tmp)
    
    ci <- tmp[,c("lower", "upper")]
  }else{
    my_coefficients <- coef(model)
    coef_names <- names(my_coefficients)
    
    ci <- suppressMessages(confint(model, level=level))
  }
  
  if (skip_intercept){
    intercept <- grep("intercept", coef_names, 
                      ignore.case = TRUE)
    if (length(intercept) > 0){
      my_coefficients <- my_coefficients[-intercept]
      ci <- ci[-intercept,]
      coef_names <- coef_names[-intercept]
    }
  }
  
  # Use the exp() if logit or cox regression
  if (inherits(model, "coxph") ||
        (!is.null(model$family$link) &&
           model$family$link %in% c("logit", "log"))){
    my_coefficients <- exp(my_coefficients)
    ci <- exp(ci)
  }
  
  if (length(my_coefficients) > 1)
    ret_val <- cbind(my_coefficients, ci)
  else
    ret_val <- matrix(c(my_coefficients, ci), nrow=1)
  
  colnames(ret_val) <- c("Estimate", 
                         sprintf("%.1f%%", 100*(1-level)/2),
                         sprintf("%.1f%%", 100*(level + (1-level)/2)))
  rownames(ret_val) <- coef_names
  return(ret_val)
}


#' @section The rms:
#' 
#' The rms-package does not have confint implemented and it is therefore a 
#' better option to go through the summary function (\code{rms:::summary.rms}).
#' Infortunately skip intercept is not an option as the summary doesn't 
#' include the intercept for the rms regression outputs
#' 
#' @param vn The variable names
#' @param data The data set
#' 
#' @importFrom Gmisc fastDoCall
#' @rdname GetCoefAndCI
prCaRmsGetCoefAndCI <- function(model, level, vn, data){
  # We need to standardize the summary call so that it
  # uses the first level as reference since this is not
  # always the case for the summary function
  # Furthermore, we need to set continuous variables to
  # 0 and 1 as this is not certain to be the case
  scall <- list(object=model,
                vnames="names",
                conf.int=level,
                est.all=FALSE)
  # Select reference for each summary call
  for(name in vn){
    if (is.factor(data[[name]])){
      scall[[name]] = levels(data[[name]])[1]
    }else if(is.logical(data[[name]])){
      # Logical should always have false as the
      # reference category
      scall[[name]] = FALSE
    }else if (length(unique(data[[name]])) == 2){
      freq <- table(data[[name]])
      scall[[name]] = sort(names(freq))[1]
      # The sort returns a string and this causes and 
      # error if the data is actually a numeric variable
      # therefore we need to change it back
      if (is.numeric(data[[name]])){
        scall[[name]] <- as.numeric(scall[[name]])
      }
    }else{
      # Perhaps a little overkill but it seems better to
      # set it to one step inside the range than just 0 vs 1
      scall[[name]] = c(median(data[[name]], na.rm=TRUE), 
                        median(data[[name]], na.rm=TRUE) + 1)
    }
  }
  # Call the antilog if this is a binomial logit
  if (!is.null(model$family$link) &&
        model$family$link == "logit")
    scall["antilog"] = TRUE
  
  s <- suppressMessages(fastDoCall(summary, scall))
  
  # Use the hazard ratios, the odds ratios or the antilog if
  # the function needs to be exponentiated
  r_no <- grepl("^ ((Hazard|Odds) Ratio|Antilog)", rownames(s))
  if (any(r_no)){
    rn <- rownames(s)[which(r_no)-1]
    s <- s[r_no, , drop=FALSE]
    rownames(s) <- rn
  }
  
  # Clean out the variable name for factors
  # rownames(s) <- gsub("^[^-]+- ([^:]+).+", "\\1", rownames(s))
  
  # Select the columns of interest
  ret_val <- s[, grep("(Effect|Lower|Upper)", colnames(s)), drop=FALSE]
  
  colnames(ret_val) <- c("Estimate", 
                         sprintf("%.1f%%", 100*(1-level)/2),
                         sprintf("%.1f%%", 100*(level + (1-level)/2)))
  return(ret_val)
}
