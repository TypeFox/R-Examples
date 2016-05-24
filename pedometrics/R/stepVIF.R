#' Variable selection using the variance-inflation factor
#' 
#' This function takes a linear model and selects the subset of predictor 
#' variables that meet a user-specific collinearity threshold measured by the 
#' variance-inflation factor (VIF).
#' 
#' @param model Linear model (object of class 'lm') containing collinear 
#' predictor variables.
#' @param threshold Positive number defining the maximum allowed VIF. Defaults 
#' to \code{threshold = 10}.
#' @param verbose Logical for indicating if iteration results should be 
#' printed. Defaults to \code{verbose = FALSE}.
#' 
#' @details
#' \code{stepVIF} starts computing the VIF of all predictor variables in the 
#' linear model. Because some predictor variables can have more than one degree 
#' of freedom, such as categorical variables, generalized variance-inflation 
#' factors (Fox and Monette, 1992) are calculated instead using 
#' \code{\link[car]{vif}}. Generalized variance-inflation factors (GVIF) 
#' consist of VIF corrected to the number of degrees of freedom (df) of the 
#' predictor variable:
#' 
#' \eqn{GVIF = VIF^{1/(2\times df)}}{GVIF = VIF^[1/(2*df)]}
#' 
#' GVIF are interpretable as the inflation in size of the confidence ellipse or 
#' ellipsoid for the coefficients of the predictor variable in comparison with 
#' what would be obtained for orthogonal data (Fox and Weisberg, 2011).
#' 
#' The next step is to evaluate if any of the predictor variables has a VIF 
#' larger than the specified threshold. Because \code{stepVIF} estimates GVIF 
#' and the threshold corresponds to a VIF value, the last is transformed to the 
#' scale of GVIF by taking its square root. If there is only one predictor 
#' variable that does not meet the VIF threshold, it is automatically removed 
#' from the model and no further processing occurs. When there are two or more 
#' predictor variables that do not meet the VIF threshold, \code{stepVIF} fits 
#' a linear model between each of them and the dependent variable. The 
#' predictor variable with the lowest adjusted coefficient of determination is 
#' dropped from the model and new coefficients are calculated, resulting in a 
#' new linear model.
#' 
#' This process lasts until all predictor variables included in the new model 
#' meet the VIF threshold.
#' 
#' Nothing is done if all predictor variables have a VIF value inferior to the 
#' threshold, and \code{stepVIF} returns the original linear model.
#' 
#' @return A linear model (object of class \sQuote{lm}) with low collinearity.
#' 
#' @references
#' Fox, J. and Monette, G. (1992) Generalized collinearity diagnostics. 
#' \emph{JASA}, \bold{87}, 178--183.
#' 
#' Fox, J. (2008) \emph{Applied Regression Analysis and Generalized Linear 
#' Models}, Second Edition. Sage.
#' 
#' Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Regression}, 
#' Second Edition. Thousand Oaks: Sage.
#' 
#' Hair, J. F., Black, B., Babin, B. and Anderson, R. E. (2010) 
#' \emph{Multivariate data analysis}. New Jersey: Pearson Prentice Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics 
#' with S.} Fourth edition. Springer.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @note The function name \code{stepVIF} is a variant of the widely used 
#' function \code{\link[MASS]{stepAIC}}.
#' 
#' @section TODO:
#' Include other criteria (RMSE, AIC, etc) as option to drop collinear 
#' predictor variables.
#' 
#' @seealso \code{\link[car]{vif}}, \code{\link[MASS]{stepAIC}}.
#' @export
#' @keywords methods regression
#' @examples
#' require(car)
#' fit <- lm(prestige ~ income + education + type, data = Duncan)
#' fit <- stepVIF(fit, threshold = 10, verbose = TRUE)
#' 
# FUNCTION #####################################################################
stepVIF <-
  function (model, threshold = 10, verbose = FALSE) {
    
    # Check if suggested packages are installed
    if (!requireNamespace("car", quietly = TRUE)) {
      stop(paste("Package 'car' needed for this function to work. ",
                 "Please install it.", sep = ""), call. = FALSE)
    }
    
    if (!inherits(model, "lm")) {
      stop ("'model' must be of class 'lm'")
    }
    if (threshold <= 0) {
      stop ("'threshold' must be a positive number")
    }
    
    # set conditional variable (number of iterations = number of variables)
    iter <- dim(stats::model.frame(model))[2]
    iter0 <- iter

        while (iter > 0) {
      iter <- iter - 1
      # calculate generalized variance-inflation factors
      init_vif <- data.frame(car::vif(model))
      var_nam <- rownames(init_vif)
      if (ncol(init_vif) == 3) {
        init_vif <- data.frame(init_vif[, 3])
        } else {
          init_vif <- data.frame(sqrt(init_vif[, 1]))
        }
      rownames(init_vif) <- var_nam
      colnames(init_vif) <- "vif"
      # get vars with vif > threshold
      lim_vif <- sqrt(threshold)
      if (max(init_vif) <= lim_vif) {
        cat("all predictor variables have a VIF lower than the threshold\n")
        break
      }
      vars_vif <- which(init_vif > lim_vif)
      df <- stats::model.frame(model)
      df_nam <- colnames(df[, -1])[vars_vif]
      if (length(df_nam) == 1) {  # only one predictor with high VIF
        # update model formula
        #out <- names(vars_vif)
        out <- df_nam
        new_form <- stats::formula(paste(". ~ .", paste(out, collapse = "-"), 
                                         sep = "-"))
        model <- stats::update(model, new_form)
        break
      }
      vars_vif <- data.frame(df[, -1][, vars_vif])
      colnames(vars_vif) <- df_nam
      # fit lm between dependent variable and every predictor variable
      new_fit <- lapply(vars_vif, function (X) stats::lm(df[, 1] ~ X))
      r2 <- unlist(sapply(new_fit, summary)["adj.r.squared", ])
      min_r2 <- which(r2 == min(r2))
      # update model formula
      out <- names(min_r2)
      new_form <- stats::formula(paste(". ~ .", paste(out, collapse = "-"), 
                                       sep = "-"))
      model <- stats::update(model, new_form)
      end_vif <- data.frame(car::vif(model))
      if (ncol(end_vif) == 3) {
        end_vif <- end_vif[, 3]
        } else {
          end_vif <- data.frame(sqrt(end_vif[, 1]))
        }
      # check VIF
      end_vif <- max(end_vif*end_vif)
      if (verbose) {  # print the results of each iteration (if requested)
        tab <- round(init_vif[names(vars_vif),]^2)
        names(tab) <- names(vars_vif)
        cat("------------------------------------------------------------\n")
        cat(paste("iteration: ", iter0 - iter, sep = ""), "\n")
        cat(paste("collinearity cases: ", length(names(tab)), sep = ""), "\n")
        cat("generalized variance inflation factor\n")
        print(tab)
        cat("correlation with dependent variable (Adjusted R squared)\n")
        print(round(r2, 4))
        cat(paste("dropped term: ", names(min_r2), "\n"))
        print(model)
      }
      if (end_vif <= threshold)  # evaluate model vif
        break
    }
    return (model)
  }
#
# End!
