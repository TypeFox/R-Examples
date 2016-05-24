print.wald.test <- function(x, ..., digits = max(3, getOption("digits") - 3)){
  Terms <- x[["Terms"]]
  b <- x[["b"]]
  H0 <- x[["H0"]]
  v <- x[["result"]][["chi2"]]
  df <- x[["df"]]
  verbose <- x[["verbose"]]
  namb <- names(b)[Terms]
  
  if(verbose){
    cat("\nCoefficients:\n")
    print(format(b, digits = digits), quote = FALSE)
    cat("\nVar-cov matrix of the coefficients:\n")
    print(format(x[["varb"]], digits = digits), quote = FALSE)
    cat("\nTest-design matrix:\n")
    print(x[["L"]])
    cat("\nPositions of tested coefficients in the vector of coefficients:",
        paste(Terms, collapse = ", "), "\n")
    if(is.null(namb))
      cat("\nH0: ",
          paste(paste(format(b[Terms], digits),
                      format(H0, digits = digits),
                      sep = " = "),
                collapse = "; "), "\n")
    else{
      cat("\nH0: ",
          paste(paste(namb, format(H0, digits = digits), sep = " = "),
                collapse = "; "), "\n")
    }
	#cat("\nTest results:\n")
  }
  if(is.null(df)) {
  	cat("\nChi-squared test:\n")
  	cat("X2 = ", format(v["chi2"], digits = digits), ", df = ", v["df"],
  		", P(> X2) = ", format(v["P"], digits = digits), "\n", sep = "")
    }
  else{
    v <- x[["result"]][["Ftest"]]
    cat("\nF test:\n")
    cat("F = ", format(v["Fstat"], digits = digits), 
        ", df1 = ", v["df1"],
        ", df2 = ", v["df2"],
        ", P(> F) = ", format(v["P"], digits = digits), "\n\n", sep = "")
    }
  }
