wald.test <- function(Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL, verbose = FALSE){
## cross checks between Terms and L
  if(is.null(Terms) & is.null(L))
    stop("One of the arguments Terms or L must be used.")
  if(!is.null(Terms) & !is.null(L))
    stop("Only one of the arguments Terms or L must be used.")
## Terms (coef of which) to be tested
  if(is.null(Terms)){
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
    }
  else
    w <- length(Terms)
## null hypothesis
  if(is.null(H0))
    H0 <- rep(0, w)
  if(w != length(H0))
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
## L matrix
  if(is.null(L)){
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for(i in 1:w)
      L[i, Terms[i]] <- 1
    }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))), sep = ""), names(b))
## computations
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if(is.null(df))
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else{
    fstat <- stat / nrow(L)
    df1 <- nrow(L); df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p),
                Ftest = c(Fstat = fstat, df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
    }
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0, L = L, result = res, verbose = verbose, df = df),
            class = "wald.test")
  }

print.wald.test <- function(x, digits = 2, ...){
  Terms <- x[["Terms"]]; b <- x[["b"]]; H0 <- x[["H0"]]; v <- x[["result"]][["chi2"]]; df <- x[["df"]]
  verbose <- x[["verbose"]]
  namb <- names(b)[Terms]
  cat("Wald test:\n", "----------\n", sep = "")
  if(verbose){
    cat("\nCoefficients:\n")
    print(format(b, digits = digits), quote = FALSE)
    cat("\nVar-cov matrix of the coefficients:\n")
    print(format(x[["Sigma"]], digits = digits), quote = FALSE)
    cat("\nTest-design matrix:\n")
    print(x[["L"]])
    cat("\nPositions of tested coefficients in the vector of coefficients:", paste(Terms, collapse = ", "), "\n")
    if(is.null(namb))
      cat("\nH0: ", paste(paste(format(b[Terms], digits), format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
    else{
      cat("\nH0: ", paste(paste(namb, format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
      }
#    cat("\nTest results:\n")
    }
  cat("\nChi-squared test:\n")
  cat("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), ", df = ", v["df"],
      ", P(> X2) = ", format(v["P"], digits = digits, nsmall = 1), "\n", sep = "")
  if(!is.null(df)){
    v <- x[["result"]][["Ftest"]]
    cat("\nF test:\n")
    cat("W = ", format(v["Fstat"], digits = digits, nsmall = 1), 
        ", df1 = ", v["df1"],
        ", df2 = ", v["df2"],
        ", P(> W) = ", format(v["P"], digits = digits), "\n", sep = "")
    }
  }
