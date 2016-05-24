K.factor <- function (n, f = NULL, alpha = 0.05, P = 0.99, side = 1, method = c("HE", 
    "HE2", "WBE", "ELL", "KM", "EXACT", "OCT"), m = 50) 
{
    if (is.null(f)) 
        f <- n - 1
    if ((length(n) != length(f)) & length(f) > 1) {
        stop(paste("Length of 'f' needs to match length of 'n'!", 
            "\n"))
    }
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    method <- match.arg(method)
    if (side == 1) {
        z.p <- qnorm(P)
        ncp <- sqrt(n) * z.p
        t.a <- suppressWarnings(qt(1 - alpha, df = f, ncp = ncp))
        K <- t.a/sqrt(n)
    }
    else {
        K.temp <- function(n, f, alpha, P, method = c("HE", "HE2", "WBE", 
            "ELL", "KM", "EXACT", "OCT"), m) {
            method <- match.arg(method)
            chi.a <- qchisq(alpha, f)
            k2 <- sqrt(f * qchisq(P, 1, 1/n)/chi.a)
            if (method == "HE") {
                TEMP4 <- function(n, f, P, alpha) {
                  chi.a <- qchisq(alpha, f)
                  z.p <- qnorm((1 + P)/2)
                  z.a <- qnorm((2 - alpha)/2)
                  df.cut <- n^2 * (1 + 1/z.a^2)
                  V <- 1 + z.a^2/n + ((3 - z.p^2) * z.a^4)/(6 * 
                    n^2)
                  K.1 <- suppressWarnings(z.p * sqrt(V * (1 + 
                    (n * V/(2 * f)) * (1 + 1/z.a^2))))
                  G <- (f - 2 - chi.a)/(2 * (n + 1)^2)
                  K.2 <- suppressWarnings(z.p * sqrt(((f * (1 + 
                    1/n))/(chi.a)) * (1 + G)))
                  if (f > df.cut) {
                    K <- K.1
                  }
                  else {
                    K <- K.2
                    if (is.na(K)) 
                      K <- 0
                  }
                  K
                }
                TEMP5 = Vectorize(TEMP4)
                K <- TEMP5(n, f, P, alpha)
            }
            else if (method == "HE2") {
                z.p <- qnorm((1 + P)/2)
                K <- z.p * sqrt((1 + 1/n) * f/chi.a)
            }
            else if (method == "WBE") {
                r <- 0.5
                delta <- 1
                while (abs(delta) > 1e-08) {
                  P.new <- pnorm(1/sqrt(n) + r) - pnorm(1/sqrt(n) - 
                    r)
                  delta <- P.new - P
                  diff <- dnorm(1/sqrt(n) + r) + dnorm(1/sqrt(n) - 
                    r)
                  r <- r - delta/diff
                }
                K <- r * sqrt(f/chi.a)
            }
            else if (method == "EXACT") {
                fun1 <- function(z, df1, P, X, n) pchisq(df1 * 
                  qchisq(P, 1, z^2)/X^2, df = df1, lower.tail = FALSE) * 
                  exp(-0.5 * n * z^2)
                fun2 <- function(X, df1, P, n, alpha, m) integrate(fun1, 
                  lower = 0, upper = 5, df1 = df1, P = P, X = X, 
                  n = n, subdivisions = m)$value
                fun3 <- function(X, df1, P, n, alpha, m) sqrt(2 * 
                  n/pi) * suppressWarnings(fun2(X, df1, P, n, 
                  alpha, m)) - (1 - alpha)
                K <- uniroot(f = fun3, interval = c(0, k2 + 1000/n), 
                  df1 = f, P = P, n = n, alpha = alpha, m = m, 
                  tol = .Machine$double.eps^0.5)$root
            }
            else if (method == "ELL") {
                if (f < (n^2)) 
                  warning("The Ellison method should only be used for f appreciably larger than n^2.", 
                    call. = FALSE)
                r <- 0.5
                delta <- 1
                z.p <- qnorm((1 + P)/2)
                while (abs(delta) > 1e-08) {
                  P.new <- pnorm(z.p/sqrt(n) + r) - pnorm(z.p/sqrt(n) - 
                    r)
                  delta <- P.new - P
                  diff <- dnorm(z.p/sqrt(n) + r) + dnorm(z.p/sqrt(n) - 
                    r)
                  r <- r - delta/diff
                }
                K <- r * sqrt(f/chi.a)
            }
            else if (method == "KM") {
                K <- k2
            }
            else if (method == "OCT") {
                delta <- sqrt(n) * qnorm((1 + P)/2)
                Fun1 <- function(z, P, ke, n, f1, delta) (2 * pnorm(-delta + 
                  (ke * sqrt(n * z))/(sqrt(f1))) - 1) * dchisq(z, f1)
                Fun2 <- function(ke, P, n, f1, alpha, m, delta) integrate(Fun1, 
                  lower = f1 * delta^2/(ke^2 * n), upper = 1000 * 
                    n, P = P, ke = ke, n = n, f1 = f1, delta = delta, 
                  subdivisions = m)$value
                Fun3 <- function(ke, P, n, f1, alpha, m, delta) abs((Fun2(ke = ke, 
                	P = P, n = n, f1 = f1, alpha = alpha, m = m, delta = delta)) - 
                  	(1 - alpha))
                K <- optim(par = k2,   fn = Fun3, lower=0, 
                  P = P, n = n, f1 = f, alpha = alpha, m = m, delta = delta, 
                  method="L-BFGS-B")$par
            }
        }
        TEMP <- Vectorize(K.temp)
        K <- TEMP(n = n, f = f, alpha = alpha, P = P, method = method, 
            m = m)
    }
    K
}
