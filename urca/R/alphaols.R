##
## alphaols
##
alphaols <- function(z, reg.number = NULL) 
{
  if (!(class(z) == "ca.jo")) {
    stop("\nPlease, provide object of class 'ca.jo' as 'z'.\n")
  }
  RKV <- z@RK%*%z@V
  colnames(RKV) <- paste("V", colnames(z@RK), sep=".")
  P <- z@P
  data.mat <- data.frame(z@R0, RKV)
  text <- colnames(data.mat)[-c(1:P)]
  text1 <- paste(text, "", sep = "+", collapse = "")
  text2 <- paste("~", substr(text1, 1, nchar(text1) - 1))
  if (!is.null(reg.number)) {
    reg.number <- as.integer(reg.number)
    if (reg.number > ncol(z@R0) || reg.number < 1) {
      stop("\nPlease, provide a valid number of the regression within \n the VECM, numbering from 1st to last row.\n")
    }
    form1 <- formula(paste("z@R0[, reg.number]", text2, "-1"))
    return(lm(substitute(form1), data = data.mat))
  }
  else if (is.null(reg.number)) {
    form1 <- formula(paste("z@R0", text2, "-1"))
    return(lm(substitute(form1), data = data.mat))
  }
}
