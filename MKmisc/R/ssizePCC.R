## gamma: tolerance between PCC(infty) and PCC(n)
## stdFC: standardized fold-change
## prev: pervalence
## nrFeatures: total number of features
## sigFeatures: number of significant features
## verbose: print 
ssize.pcc <- function(gamma, stdFC, prev = 0.5, nrFeatures, sigFeatures = 20, verbose = FALSE){
  pcc.n.a <- function(alpha, n, d.s, nrFeatures, m){
    pow <- pt(d.s*sqrt(n)+qt(alpha/2, df = n-2), df = n-2)
    T1 <- d.s*m*pow
    T2 <- sqrt(m*pow + alpha*(nrFeatures-m))
    pnorm(T1/T2)
  }
  d.s <- stdFC/2
  n <- 2
  M <- sigFeatures
  repeat{
    n <- n + 1
    m <- seq_len(M)
    alpha <- numeric(M)
    for(mi in m){
      alpha[mi] <- optimize(f = pcc.n.a, interval = c(0, 1), maximum = TRUE,
                            n = n, d.s = d.s, nrFeatures = nrFeatures, m = mi,
                            tol = 1e-10)$maximum
    }
    pow <- pt(d.s*sqrt(n)+qt(alpha/2, df = n-2), df = n-2)
    T1 <- d.s*m*pow
    T2 <- sqrt(m*pow + alpha*(nrFeatures-m))
    crit <- max(pnorm(d.s*sqrt(m)) - pnorm(T1/T2))
    if(verbose){
      cat("========================================\n")
      cat("Current sample size:\t", n, "\n")
      cat("Upper bound Un:\t", crit, "\n")
    }
    if(crit < gamma) break
  }
  n1 <- 0.5*n/min(prev, 1-prev)
  ns <- ceiling(c(prev*n1, (1-prev)*n1))
  NOTE <- "n1 is number of cases, n2 is number of controls"
  METHOD <- "Sample Size Planning for Developing Classifiers Using High Dimensional Data"
  
  res <- structure(list(gamma = gamma, prev = prev,
                        nrFeatures = nrFeatures,
                        n1 = ns[1], n2 = ns[2],
                        note = NOTE, method = METHOD), 
                   class = "power.htest")  
  res
}
