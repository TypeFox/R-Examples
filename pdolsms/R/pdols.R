#' pdolsms: Cointegration Vector Estimation by Panel DOLS
#'
#' The pdolsms package provides a function for the estimation of panel data cointegrating relationships following Mark and Sul (2003), <DOI:10.1111/j.1468-0084.2003.00066.x>.
#'
#' @docType package
#' @name pdolsms
NULL
#> NULL


#' Cointegration Vector Estimation by Panel DOLS
#'
#' \code{pdols} estimates panel data cointegrating relationships following the estimator of Mark and Sul (2003), <DOI:10.1111/j.1468-0084.2003.00066.x>.
#' @param formula An object of class formula. Note, pdols does not know how to deal with factors.
#' @param data A data frame in long panel format where each period of data is a seperate row for each individual.
#' @param index A vector of two strings for the individual and the time column names.
#' @param p An integer indicating the number of leads and lags to include in the regression.
#' @param icase One of the following character strings "noconstant", "constant", or "trend". Specifying "trend" also includes a constant in the estimation.
#'
#' @return \code{pdols} returns an object of class pdols with methods print, summary and coef.
#' This object contains a list with the dynamic OLS coefficients for each individual, as well as pooled ordinary and common time effects dols estimates of the coefficients. For more details, please refer to Mark and Sul (2003), <DOI:10.1111/j.1468-0084.2003.00066.x>.
#' Two standard errors are reported. The first is based on Andrews and Monahan's Pre-whitening methods and the second is based on parametric correction.
#' 
#' @references Mark, N. C. and Sul, D. (2003), Cointegration Vector Estimation by Panel DOLS and Long-run Money Demand. Oxford Bulletin of Economics and Statistics, 65: 655-680. <DOI:10.1111/j.1468-0084.2003.00066.x>
#' @export
#' @import stats
#' @importFrom reshape2 dcast
#' @examples
#' data("marksul2003")
#' result <- pdols(formula = Y ~ X + Z,
#'                 index = c("Country", "Year"),
#'                 data = marksul2003, 
#'                 p = 2, 
#'                 icase = "constant")
#' print(result)                 

pdols <- function(formula, data, index, p, icase = c("noconstant", "constant", "trend")) {
  
  cl <- match.call()
  args <- list(formula = formula, index = index, p = p, icase = icase)
  
  icase <- match.arg(icase)
  icase <- switch(icase,
                  noconstant = 1,
                  constant = 2,
                  trend = 3)
  
  
  if (class(data) != "data.frame") stop("Data must be specified as a data.frame")
  if (!all(all.vars(formula) %in% names(data))) stop("Variables specified in the formula are not found in the dataset")
  if (any(sapply(data[, all.vars(formula) ], is.factor))) stop("dols does not know how to deal with factors")
  if (p < 0) stop("Number of lags and leads, p, must be greater than or equal to 0")
  if (!is.character(index) | (length(index) < 2)) stop("Please specify a vector of two strings for the individual and the time column names")
  if (!all(table(data[,index]) == 1)) stop("The data supplied must be a balanced panel with one observation for each individual for each period")
  if (any(is.na(data))) stop("The data supplied must be a balanced panel with one observation for each individual for each period, pdols does not know how to deal with missing values.")
  if (!is.numeric(data[, index[2]])) stop("The time variable must be numeric. Please coerce it to this class.")
  
  data <- cbind(data[, c(index, all.vars(formula)[1])], model.matrix(formula, data))
  colnames(data)[1:2] <- index
  data <- data[, !colnames(data) %in% "(Intercept)"]
  
  data <- data[order(data[index[2]], data[index[1]]), ]
  
  ldata <- vector(mode = "list", length = NCOL(data) - 2)
  names(ldata) <- colnames(data)[3:NCOL(data)]  
  
  cn <- colnames(data)
  for (i in 1:length(ldata)) {
    vn <- cn[i + 2]
    cformula <- paste(cn[2], "~", cn[1])
    
    ldata[[i]] <- as.matrix(dcast(data, formula(cformula), value.var = vn))
    rownames(ldata[[i]]) <- ldata[[i]][, 1]
    ldata[[i]] <- ldata[[i]][, -1]
  }
  
  if (NROW(ldata[[1]]) < 2*p + 1) stop("Number of periods is less than the ottal number of lags and leads")
  if (NCOL(ldata[[1]]) < 2) stop("There are less than two individuals in the panel")
  
  
  results <- pdols.fit(ldata, p = p, icase = icase)
  
  
  # Need to do some renaming for prettier printing
  rownames(results$`Pooled OLS ORD`) <- colnames(data)[4:NCOL(data)]
  rownames(results$`Pooled OLS CTE`) <- colnames(data)[4:NCOL(data)]
  
  colnames(results$`DOLS Single Equation`) <- colnames(ldata[[1]])
  rownames(results$`DOLS Single Equation`)[rownames(results$`DOLS Single Equation`) 
                                           %in% paste0("A", 1:(NCOL(data)-2))] <- colnames(data)[4:NCOL(data)]
  if (NCOL(data) == 4) {
    results$`DOLS Single Equation` <- results$`DOLS Single Equation`[-3, ]
  }
  
  
  results <- append(results, list(call = cl, args = args))
  class(results) <- "pdols"
  return(results)
}






#' Panel DOLS Coefficients
#'
#' @param object An object of class pdols
#' @param type The type of coefficient to return. Either ordinary dols ("ordinary") or common time effects ("cte")
#'
#' @return A vector of the specified coefficients
#' @param ...	Further arguments passed to or from other methods.
#' @export
#' @import stats
#' @references Mark, N. C. and Sul, D. (2003), Cointegration Vector Estimation by Panel DOLS and Long-run Money Demand. Oxford Bulletin of Economics and Statistics, 65: 655-680. <DOI:10.1111/j.1468-0084.2003.00066.x>
#' @examples
#' data("marksul2003")
#' result <- pdols(formula = Y ~ X + Z - 1,
#'                 index = c("Country", "Year"),
#'                 data = marksul2003, 
#'                 p = 2, 
#'                 icase = "constant")
#' coef(result, type = "cte") 
#' coef(result, type = "o") 


coef.pdols <- function(object, type = c("ordinary", "cte"), ...) {
  
  type <- match.arg(type)
  coef <- switch(type,
                 ordinary = setNames(object$`Pooled OLS ORD`[, 1], 
                                     rownames(object$`Pooled OLS ORD`)),
                 cte = setNames(object$`Pooled OLS CTE`[, 1],
                                rownames(object$`Pooled OLS CTE`))
  )
  
  if (length(coef)) {
    return (coef)
  } else {
    stop("Coefficients not found")
  }
}




#' A print method for the pdols class
#'
#' @param x An object of class pdols
#' @param type The type of coefficient to return. Either ordinary dols ("ordinar") or common time effects ("cte")
#' @param digits The number of digits to display when printing the coefficients.
#' @param ...	Further arguments passed to or from other methods.
#' @export
#' @import stats
#' @references Mark, N. C. and Sul, D. (2003), Cointegration Vector Estimation by Panel DOLS and Long-run Money Demand. Oxford Bulletin of Economics and Statistics, 65: 655-680. <DOI:10.1111/j.1468-0084.2003.00066.x>
#' @examples
#' data("marksul2003")
#' result <- pdols(formula = Y ~ X + Z - 1,
#'                 index = c("Country", "Year"),
#'                 data = marksul2003, 
#'                 p = 2, 
#'                 icase = "constant")
#' print(result) 

print.pdols <- function(x, type = c("ordinary", "cte"), digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call, quote = FALSE)
  cat("\n")
  
  
  type <- match.arg(type)
  s <- switch(type,
              ordinary = "Ordinary DOLS Coefficients:",
              cte = "Common Time Effects Coefficients:"
  )
  
  if (length(coef(x, type))) {
    cat(s, "\n")
    print(format(coef(x, type), digits = digits), print.gap = 2L, quote = FALSE)
  } else  {
    cat("Coefficients not found")
  }
  cat("\n\n") 
  invisible(x)
}




#' A summary method for the pdols class
#'
#' @param object An object of class pdols
#' @param ...	Further arguments passed to or from other methods.
#' @export
#' @import stats
#' @references Mark, N. C. and Sul, D. (2003), Cointegration Vector Estimation by Panel DOLS and Long-run Money Demand. Oxford Bulletin of Economics and Statistics, 65: 655-680. <DOI:10.1111/j.1468-0084.2003.00066.x>
#' @examples
#' data("marksul2003")
#' result <- pdols(formula = Y ~ X + Z - 1,
#'                 index = c("Country", "Year"),
#'                 data = marksul2003, 
#'                 p = 2, 
#'                 icase = "constant")
#' summary(result) 

summary.pdols <- function(object, ...) {
  
  icase <- match.arg(object$args$icase, choices = c("noconstant", "constant", "trend"))
  icase <- switch(icase,
                  noconstant = 1,
                  constant = 2,
                  trend = 3)
  
  
  cat("
      ________________________________________________
      Reference: 
      Mark, Nelson C., and Donggyu Sul. 
      'Cointegration Vector Estimation by 
      Panel DOLS and Long-run Money Demand.' 
      Oxford Bulletin of Economics and Statistics 65.5 
      (2003): 655-680.
      ________________________________________________\n\n")
  
  cat("S.E. Andrews : Standard error based on Andrews and Monahan's Pre-whitening method\n")
  cat("S.E. Param   : Standard error based on parametric correction\n\n\n")
  
  cat("Call:\n")
  print(object$call, quote = FALSE)
  cat("\n\n\n")
  
  if (icase == 2) {
    cat("Model: y(it) = c(i) + A*X(it) + u(it)\n\n\n")
  }
  if (icase == 3) {
    cat("Model: y(it) = c(i) + d(i)t + A*X(it) + u(it)\n\n\n")
  }
  
  cat("Dynamic OLS for Single Equation\n\n")
  print(object$`DOLS Single Equation`)
  cat("\n\n")
  
  cat("Pooled OLS - Common Time Effects\n\n")
  print(object$`Pooled OLS CTE`)
  cat("\n\n")
  
  cat("Pooled OLS - Ordinary\n\n")
  print(object$`Pooled OLS ORD`)
  cat("\n\n")
  
  invisible(object)
}

#' Sample panel data for pdols
#'
#' Sample panel data in long format for pdols
#'
#' @format A data frame with 760 rows and 5 variables:
#' @source \url{http://www.utdallas.edu/~d.sul/papers/}
"marksul2003"
