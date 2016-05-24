wald.test <- function(b, varb, Terms = NULL, L = NULL, H0 = NULL, df = NULL, verbose = FALSE, ...){
  
	# cross checks between Terms and L
  if(is.null(Terms) & is.null(L))
    stop("One of the arguments Terms or L must be used.")
  if(!is.null(Terms) & !is.null(L))
    stop("Only one of the arguments Terms or L must be used.")

	# Terms (coef of which) to be tested
  if(is.null(Terms)){
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
    }
  else
    w <- length(Terms)

	# null hypothesis
  if(is.null(H0))
    H0 <- rep(0, w)
  if(w != length(H0))
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")

	# L matrix
  if(is.null(L)){
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for(i in 1:w)
      L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))), sep = ""), names(b))

	# computations
  f <- L %*% b
  V <- varb
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

	structure(list(varb = varb, b = b, Terms = Terms, H0 = H0, L = L,
                 result = res, verbose = verbose, df = df),
            class = "wald.test")

}
