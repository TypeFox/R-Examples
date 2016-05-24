context("Consistency between sparse/non-sparse encoding")

test_that("Consistency of quadrupen between sparse/non-sparse encoding of the predictors", {

  ## data generation
  rlm.sparse <- function(n, p, reg="low", mu=3, size=5, prob=0.01) {
    s <- switch(reg,
                low  = floor(0.50 * min(n,p)),
                med  = floor(0.10 * min(n,p)),
                high = floor(0.01 * min(n,p)))

    w <- rep(0,p)
    w[sample(1:p,s)] <- sample(c(-1,1),s,replace=TRUE)*runif(s,1,2)

    X <- matrix(rbinom(n*p, size, prob),n,p)

    Xw <- crossprod(t(X),w)
    sigma <- sum(Xw^2) * 0.01/n
    epsilon <- rnorm(n) * sigma
    y <- mu + Xw + epsilon
    r2 <- 1 - sum(epsilon^2) / sum((y-mean(y))^2)

    return(list(y=y, x=X, r2=r2, w=w, mu=mu))
  }

  n <- 500
  p <- 10000
  lambda2 <- 0.05
  data <- rlm.sparse(n, p, prob=0.01, reg="med")
  s <- sum(data$w !=0)
  y    <- data$y
  x.ns <- as.matrix(data$x)
  x.sp <- Matrix(data$x, sparse=TRUE)
  max.feat <- s * 10

  cat("\n\tProblem with",p, "predictors and",n,"samples,",s,"true nonzeros in beta.star.")
  cat("\n\tThe densely encoded design matrix weights" , round(object.size(x.ns) /(1024^2),2), "Mo.")
  cat("\n\tThe sparsely encoded design matrix weights", round(object.size(x.sp) /(1024^2),2), "Mo.")

  cat("\n\tdense coding...")
  out.enet.ns <- elastic.net(x.ns, y, lambda2=lambda2, max.feat=max.feat, min.ratio=1e-3, control=list(timer=TRUE))

  cat(" took", out.enet.ns@monitoring$external.timer, "seconds to activate",
      rowSums(out.enet.ns@active.set)[length(out.enet.ns@lambda1)],"variables.")

  cat("\n\tsparse coding...")
  out.enet.sp <- elastic.net(x.sp, y, lambda2=lambda2, max.feat=max.feat, min.ratio=1e-3, control=list(timer=TRUE))
  cat(" took", out.enet.sp@monitoring$external.timer, "seconds to activate",
      rowSums(out.enet.sp@active.set)[length(out.enet.sp@lambda1)],"variables.")

  ## remove timer for fair comparison!!!
  out.enet.sp@monitoring[5:7] <- NULL
  out.enet.ns@monitoring[5:7] <- NULL
  expect_that(out.enet.sp, equals(out.enet.ns))
})


