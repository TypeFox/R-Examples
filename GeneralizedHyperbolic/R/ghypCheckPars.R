### Check for boundary and error cases
### Return case normal, or error
### In case of error, include appropriate message
ghypCheckPars <- function(param) {

  param <- as.numeric(param)
  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  case <- ""
  errMessage <- ""

  if (length(param) != 5) {
    case <- "error"
    errMessage <- "param vector must contain 5 values"
  } else {
    if (alpha < 0) {
      case <- "error"
      errMessage <- "alpha must not be less than zero"
    } else {
      if (lambda == 0) {
        if (abs(beta) >= alpha) {
          case <- "error"
          errMessage <-
              "absolute value of beta must be less than alpha when lambda = 0"
        }

        if (delta <= 0) {
            case <- "error"
            errMessage <- "delta must be greater than zero when lambda = 0"
        }

        if (abs(beta) >= alpha & delta <= 0) {
            case <- "error"
            errMessage <- "absolute value of beta must be less than alpha and delta must be greater than zero when lambda = 0"
        }
      }

      if (lambda > 0) {
        if (abs(beta) >= alpha) {
          case <- "error"
          errMessage <- "absolute value of beta must be less than alpha when lambda > 0"
        }

        if (delta < 0) {
          case <- "error"
          errMessage <- "delta must be non-negative when lambda > 0"
        }

        if (abs(beta) >= alpha & delta < 0) {
          case <- "error"
          errMessage <- "absolute value of beta must be less than alpha and delta must be less than zero when lambda > 0"
        }

        if (case != "error") {
          if (lambda == 1) {
            if (alpha > 0 & abs(beta) < abs(alpha) & delta == 0)
              case <- "skew laplace"

            if (alpha > 0 & beta == 0 & delta == 0)
              case <- "laplace"

            if (alpha > 0 & abs(beta) < abs(alpha) & delta > 0)
              case <- "hyperbolic"
          } else {
            if (alpha > 0 & abs(beta) < abs(alpha) & delta == 0)
              case <- "variance gamma"
          }
        }
      }

      if (lambda < 0) {
        if (abs(beta) > alpha) {
          case <- "error"
          errMessage <- "absolute value of beta must be less than or equal to alpha when lambda < 0"
        }

        if (delta <= 0) {
          case <- "error"
          errMessage <- "delta must be greater than zero when lambda < 0"
        }

        if (abs(beta) > alpha & delta <= 0) {
          case <- "error"
          errMessage <-
              "absolute value of beta must be less than or equal to alpha and delta must be greater than zero when lambda < 0"
        }

        if (case != "error") {
          if (lambda == -1/2) {
            if (alpha == 0 & beta == 0 & delta > 0)
              case <- "cauchy"

            if (alpha > 0 & abs(beta) < abs(alpha) & delta > 0)
              case <- "normal inverse gaussian"
          } else {
            # Allowing a tolerance for the alpha == abs(beta) test
            if (abs(alpha - abs(beta)) < 0.001 & beta >= 0 & delta > 0)
              case <- "skew hyperbolic"

            if (alpha == 0 & beta == 0 & delta > 0)
              case <- "student's t"
          }
        }
      }
    }
  }

  # Using the normal distribution as a limiting case when no distribution
  # has been set and no errors have been found
  if (case == "")
    case <- "normal"

  result <- list(case = case, errMessage = errMessage)
  return(result)
}
