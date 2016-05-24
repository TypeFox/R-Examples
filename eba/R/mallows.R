mbt <- function(data, bootstrap = FALSE, nsim = 1000, ...){
  # Estimates parameters of Mallows-Bradley-Terry model
  # data: data.frame with columns 1:t being the ranks, column t + 1 the counts
  # author: Florian Wickelmaier (wickelmaier@web.de)
  # Last mod: Aug/31/2011
  #           Mar/27/2012 now ordering of data is irrelevant
  #                       add parametric bootstrap

  t <- ncol(data) - 1
  P <- permutations(t, t)  # need all factorial(t) permutations in X
  X <- t - P
  y <- rep(0, nrow(X))
  names(y)          <- apply(P,          1, paste, collapse = ":")
  rownames(data)    <- apply(data[,1:t], 1, paste, collapse = ":")
  y[rownames(data)] <- data[,t + 1]

  # perm.idx <- which(apply(P,          1, paste, collapse=":") %in%
  #                   apply(data[,1:t], 1, paste, collapse=":"))
  # y[perm.idx] <- data[,t + 1]
  
  glm1 <- glm(y ~ X, poisson)
  u    <- exp(c(glm1$coef[2:t], 0))  # ratio scale, normalization arbitrary
  u    <- u/sum(u)
  names(u) <- colnames(data)[1:t]

  gof <- c(deviance(glm1), glm1$df.resid, 1 - pchisq(deviance(glm1),
    glm1$df.resid), NA)
  names(gof) <- c("G2", "df", "pval", "bstp pval")

  ## Parametric bootstrap
  if(bootstrap) gof["bstp pval"] <-
    mean(sapply(simulate(glm1, nsim=nsim, ...),
      function(y1) deviance(glm.fit(X, y1, family=poisson()))) > gof["G2"])

  out <- list(coefficients=u, goodness.of.fit=gof, perm.idx=rownames(data),
    y=y, mbt.glm=glm1)
  class(out) <- "mbt"
  out
}


print.mbt <- function(x, digits=max(3, getOption("digits")-3),
  na.print="", ...){
  cat("\nMallows-Bradely-Terry (MBT) models\n\n")
  cat("Parameter estimates:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
      quote = FALSE)
  G2    <- x$goodness.of.fit[1]
  df    <- x$goodness.of.fit[2]
  pval  <- x$goodness.of.fit[3]
  bpval <- x$goodness.of.fit[4]
  cat("\nGoodness of fit (-2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits),
      ", p = ", format(pval, digits=digits),
      if(is.na(bpval)) NULL else ", bootstrap p = ",
      if(is.na(bpval)) NULL else format(bpval, digits=digits), "\n", sep="")
  cat("\n")
  invisible(x)
}


## permutations() is taken from the gtools package by Gregory R. Warnes:

##
## Original version by Bill Venables and cited by Matthew
## Wiener (mcw@ln.nimh.nih.gov) in an email to R-help dated
## Tue, 14 Dec 1999 09:11:32 -0500 (EST) in response to
## Alex Ahgarin <datamanagement@email.com>
##

permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
      {
        if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
          inner  <-  Recall(n, r-1, v)
          cbind( rep( v, rep(nrow(inner),n)  ),
                 matrix( t(inner), ncol=ncol(inner), nrow=nrow(inner) * n ,
                        byrow=TRUE )
                )
        }
      }
  else
    sub <- function(n, r, v)
      {
        if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
        X  <-  NULL
        for(i in 1:n)
          X  <-  rbind( X, cbind( v[i], Recall(n-1, r - 1, v[-i])))
        X
        }
      }

  sub(n, r, v[1:n])
}
