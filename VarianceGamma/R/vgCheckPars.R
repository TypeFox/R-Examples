### Check for error cases

vgCheckPars <- function(param, ...){
    param <- as.numeric(param)
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    case <- "normal"
    errMessage <- ""
    if (length(param) != 4) {
      case <- "error"
      errMessage <- "param vector must contain 4 values"
    } else if (sigma <= 0 ) {
         case <- "error"
         errMessage <- "sigma must be greater than zero"
     } else if (nu <= 0 ) {
         case <- "error"
         errMessage <- "nu must be greater than zero"
     }

     result <- list(case = case, errMessage = errMessage)
     return(result)
}
