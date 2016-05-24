.gofCopulapb = function (copula, x, N = 1000, method = eval(formals(gofTstat)$method), estim.method = eval(formals(fitCopula)$method), ...) 
{
  stopifnot(is(copula, "copula"), N >= 1)
  if (!is.matrix(x)) 
    x <- rbind(x, deparse.level = 0L)
  stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == 
              d)
  estim.method <- match.arg(estim.method)

  if (missing(estim.method) && !missing(method)) {
    eMeth <- eval(formals()$estim.method)
    if (!is.na(i <- pmatch(method, eMeth))) {
      warning("old (pre 0.999-*) argument 'method' is now called 'estim.method'")
      estim.method <- eMeth[i]
      method <- "Sn"
    }
  }
  
  method <- match.arg(method)
  
  uhat <- x
  C.th.n <- fitCopula(copula, uhat, method = estim.method, 
                      estimate.variance = FALSE, ...)@copula

  u = uhat
  Tstat <- if (method == "Sn") 
    gofTstat(u, method = method, copula = C.th.n)
  else gofTstat(u, method = method)
  T0 <- vapply(1:N, function(k) {
    repeat {
      Uhat <- rCopula(n, C.th.n)
      C.th.n. <- try(fitCopula(copula, Uhat, method = estim.method, 
                           estimate.variance = FALSE, ...)@copula, silent = T)
      if (class(C.th.n.) != "try-error"){break}
    }
    u. = Uhat
    T0. <- if (method == "Sn") 
      gofTstat(u., method = method, copula = C.th.n.)
    else gofTstat(u., method = method)
    T0.
  }, NA_real_)
  switch(method,
         SnB = {matrix_names = "RosenblattSnB"},
         SnC = {matrix_names = "RosenblattSnC"},
         AnChisq = {matrix_names = "ADChisq"},
         AnGamma = {matrix_names = "ADGamma"},
         Sn = {matrix_names = "Sn"})
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test with %s test", 
                                                   matrix_names),
                 erg.tests = matrix(c((sum(T0 >= Tstat) + 0.5)/(N + 1), Tstat), ncol = 2, 
                                    dimnames = list(matrix_names, c("p.value", "test statistic")))))
  
}


# .gofCopulamult = function (copula, x, N = 1000, method = "Rn", estim.method = eval(formals(fitCopula)$method), m, zeta.m, b, ...) 
# {
#   stopifnot(is(copula, "copula"), N >= 1)
#   if (!is.matrix(x)) 
#     x <- rbind(x, deparse.level = 0L)
#   stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == 
#               d)
#   estim.method <- match.arg(estim.method)
#   
#   if (missing(estim.method) && !missing(method)) {
#     eMeth <- eval(formals()$estim.method)
#     if (!is.na(i <- pmatch(method, eMeth))) {
#       warning("old (pre 0.999-*) argument 'method' is now called 'estim.method'")
#       estim.method <- eMeth[i]
#       method <- "Sn"
#     }
#   }
#   stopifnot(is(copula, "copula"), N >= 1)
#   if (!is.matrix(x)) 
#     x <- rbind(x, deparse.level = 0L)
#   stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == 
#               d)
#   method <- match.arg(method)
#   estim.method <- match.arg(estim.method)
#   if (estim.method == "ml") 
#     stop("estim.method='ml' not available")
#   if (estim.method %in% c("irho", "itau") && d > 2) 
#     stop("only bivariate case possible for estim.method='irho' or ='itau'")
#   u. <- pobs(x)
#   C.th.n <- fitCopula(copula, u., method = estim.method, estimate.variance = FALSE, 
#                       ...)@copula
#   Rn = {
#     if (estim.method != "itau") stop("Currently only estim.method='itau' available for method='Rn'")
#     C.th.n. <- pCopula(u., C.th.n)
#     denom <- (C.th.n. * (1 - C.th.n.) + zeta.m)^m
#     Cn. <- C.n(u., U = u.)
#     Tstat <- sum(((Cn. - C.th.n.)/denom)^2)
#     Z <- matrix(rnorm(N * n), nrow = N, ncol = n)
#     Zbar <- rowMeans(Z)
#     factor <- (Z - Zbar)/sqrt(n)
#     ind <- vapply(1:n, function(i) colSums(t(u.) <= u.[i, 
#                                                        ]) == d, logical(n))
#     hCnh. <- factor %*% ind
#     for (j in 1:d) {
#       u.j <- matrix(1, nrow = n, ncol = d)
#       u.j[, j] <- u.[, j]
#       ind.u.j <- vapply(1:n, function(i) colSums(t(u.) <= 
#                                                    u.j[i, ]) == d, logical(n))
#       Cnh <- factor %*% ind.u.j
#       Cjn. <- dCn(u., U = u., j.ind = j, b = b)
#       Cjn.. <- matrix(rep(Cjn., each = N), nrow = N, ncol = n)
#       hCnh. <- hCnh. - Cjn.. * Cnh
#     }
#     J.th.n <- Jscore(C.th.n, u = u., method = estim.method)
#     Theta <- rowSums(Z * rep(J.th.n, each = N)/sqrt(n))
#     d.dth.C.th.n <- as.vector(dCdtheta(C.th.n, u = u.))
#     num <- t(hCnh. - outer(Theta, d.dth.C.th.n))
#     T0 <- colSums((num/denom)^2)/n
#   }
# structure(class = "htest", list(method = sprintf("Multiplier bootstrap goodness-of-fit test with 'method'=\"%s\", 'estim.method'=\"%s\"", 
#                                                  method, estim.method), parameter = c(parameter = C.th.n@parameters), 
#                                 statistic = c(statistic = Tstat), p.value = (sum(T0 >= Tstat) + 
#                                                                                0.5)/(N + 1)))#, data.name = deparse(substitute(x))))
# }