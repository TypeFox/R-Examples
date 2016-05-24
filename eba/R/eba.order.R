# 2011/FEB/28: new zap.ind, tst.ind in printCoefmat in print.summary.eba.order
#
# 2011/MAR/01: removed dependencies on $estimate and $se components in objects
#              or class eba.order

eba.order <- function(M1, M2 = NULL, A = 1:I, s = c(rep(1/J, J), 1),
  constrained = TRUE){
  # See OptiPt
  # M1, M2: paired-comparison matrices in both within-pair orders
  # author: Florian Wickelmaier (wickelmaier@web.de)
  #
  # Fix warning: take sqrt of last element of diag(solve(hes)) only 
  #
  # last mod: 01/Mar/2011

  if(is.null(M2)){  # support 3d array
    M2 <- M1[,,2]
    M1 <- M1[,,1]
  }

  I <- ncol(M1)        # number of alternatives/stimuli
  J <- max(unlist(A))  # number of eba parameters

  idx1 <- idx0 <- matrix(0, I*(I-1)/2, J)  # index matrices
  rdx  <- 1
  for(i in 1:(I - 1)){
    for(j in (i + 1):I){
      idx1[rdx, setdiff(A[[i]], A[[j]])] <- 1
      idx0[rdx, setdiff(A[[j]], A[[i]])] <- 1
      rdx <- rdx + 1
    }
  }

  y1 <- c( t(M1)[lower.tri(M1)], t(M2)[lower.tri(M2)] )  # response vectors
  y0 <- c( M1[lower.tri(M1)], M2[lower.tri(M2)] )
  n  <- y1 + y0
  names(y1) <- names(y0) <- names(n) <- NULL
  logL.sat <- sum(dbinom(y1, n, y1/n, log=TRUE))  # likelihood of sat. model

  if(constrained){  # minimization
    out <- nlm(L.constrained.order, s, y1=y1, m=n, i1=idx1, i0=idx0)  # constr
  }else{
    out <- nlm(L.order, s, y1=y1, m=n, i1=idx1, i0=idx0)  # unconstrained
  }

  p <- out$estimate  # optimized parameters
  names(p) <- c(1:J, 'order')
  hes <- nlme::fdHess(p, L.order, y1, n, idx1, idx0)$H  # numerical Hessian
  ## Discard the order effect here
  # cova <- solve(rbind(cbind(hes, 1), c(rep(1, J+1),0)))[1:(J+1),1:(J+1)]
  cova <- solve(rbind(cbind(hes[-(J+1),-(J+1)], 1), c(rep(1, J), 0)))[1:J,1:J]
  ## Add the SE of the order effect to the se vector
  # se <- sqrt(diag(cova))  # standard error
  # se <- c(sqrt(diag(cova)), sqrt( diag(solve(hes))[J + 1] )) # standard error
  # ci <- qnorm(.975)*se    # 95% confidence interval
  ## Make block diagonal cov matrix (add zearos)
  cova <- rbind(cbind(cova, 0), c(rep(0, J), diag(solve(hes))[J + 1]))
  dimnames(cova) <- list(names(p), names(p))
  logL.eba <- -out$min    # likelihood of the specified model

  fitted1 <- fitted2 <- matrix(0, I, I)  # fitted PCMs
  eba.p <- p[-length(p)]
  n1 <- n[1:(length(n)/2)]
  n2 <- n[(length(n)/2 + 1):length(n)]
  fitted1[lower.tri(fitted1)] <- n1/(1+p["order"]*idx0%*%eba.p/idx1%*%eba.p)
  fitted2[lower.tri(fitted2)] <- n2/(1+idx0%*%eba.p/(p["order"]*idx1%*%eba.p))
  fitted1 <- t(fitted1); fitted2 = t(fitted2)
  fitted1[lower.tri(fitted1)] <- n1/(1+idx1%*%eba.p/(p["order"]*idx0%*%eba.p))
  fitted2[lower.tri(fitted2)] <- n2/(1+p["order"]*idx1%*%eba.p/idx0%*%eba.p)
  fitted <- array(c(fitted1, fitted2), c(nrow(fitted1), ncol(fitted1), 2))
  dimnames(fitted) <- list(">"=rownames(M1), "<"=colnames(M1),
    order=c("1","2"))

  # predicted probabilities
  mu <- as.numeric( c( 1 / (1 + p["order"]*idx0%*%eba.p / idx1%*%eba.p),
                       1 / (1 + idx0%*%eba.p / (p["order"]*idx1%*%eba.p)) ))

  G2   <- 2*(logL.sat - logL.eba)  # G2 goodness-of-fit statistic
  df   <- I*(I - 1) - (J + 1 - 1)
  pval <- 1 - pchisq(G2, df)
  gof  <- c(G2, df, pval)
  names(gof) <- c("-2logL", "df", "pval")
  X2   <- sum((M1 - fitted1)^2/fitted1, (M2 - fitted2)^2/fitted2, na.rm=TRUE)

  u <- numeric()  # scale values
  for(i in seq_len(I)) u = c(u, sum(p[A[[i]]]))
  names(u) <- colnames(M1)

  z <- list(coefficients=p, estimate=p, fitted=fitted,
           logL.eba=logL.eba, logL.sat=logL.sat, goodness.of.fit=gof,
           u.scale=u, hessian=-hes, cov.p=cova, idx1=idx1, idx0=idx0,
           chi.alt=X2, A=A, y1=y1, y0=y0, n=n, mu=mu, M1=M1, M2=M2)
  class(z) <- "eba.order"
  z
}


L.order <- function(p, y1, m, i1, i0)
  -sum(dbinom( y1, m,
     c( 1 / (1 + p[length(p)] * i0%*%p[-length(p)] / i1%*%p[-length(p)]),
        1 / (1 + i0%*%p[-length(p)] / (p[length(p)] * i1%*%p[-length(p)]))
      ), log=TRUE ))


L.constrained.order <- function(p, y1, m, i1, i0)  # constrain search space
  ifelse(all(p > 0),
   -sum(dbinom( y1, m,
     c( 1 / (1 + p[length(p)] * i0%*%p[-length(p)] / i1%*%p[-length(p)]),
        1 / (1 + i0%*%p[-length(p)] / (p[length(p)] * i1%*%p[-length(p)]))
      ), log=TRUE )), 1e20)


summary.eba.order <- function(object, ...){
  # Last mod: Mar/01/2011, FW

  x <- object
  I <- length(x$A)         # number of stimuli
  npairs <- I*(I - 1)/2    # number of pairs
  J <- max(unlist(x$A))    # number of eba parameters
  y <- c(x$y1, x$y0)       # response vector
  y1.1 <- x$y1[1:npairs]
  n1 <- x$n[1:npairs]
  y1.2 <- x$y1[(npairs + 1):(2*npairs)]
  n2 <- x$n[(npairs + 1):(2*npairs)]

  # coef <- x$estimate[1:J]  # only eba parameters enter in table
  coef <- coef(x)[1:J]  # only eba parameters enter in table
  # s.err <- x$se[1:J]
  s.err <- sqrt(diag(vcov(x)))[1:J]
  tvalue <- coef / s.err
  pvalue <- 2 * pnorm(-abs(tvalue))
  coef.table <- cbind(coef, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef), c("Estimate", "Std. Error",
                                              "z value", "Pr(>|z|)"))

  K <- sum(x$y0[1:npairs] + x$y1[(npairs + 1):(2*npairs)])
  N <- sum(x$n)
  order0 <- K / (N - K)  # order effect if all treatment effects are equal
  se0 <- sqrt( -N / (K * (K - N)) )  # standard error of order0

  # order.fx <- c(x$estimate["order"], order0)
  order.fx <- c(coef(x)["order"], order0)
  # order.se <- c(x$se[length(x$se)], se0)
  order.se <- c(sqrt(diag(vcov(x)))[J + 1], se0)
  o.tvalue <- (order.fx - 1) / order.se  # test against 1
  o.pvalue <- 2 * pnorm(-abs(o.tvalue))
  order.table <- cbind(order.fx, order.se, o.tvalue, o.pvalue)
  dimnames(order.table) <- list(c("order", "order0"),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  eba.p <- OptiPt(x$M1 + x$M2, x$A)  # EBA for pooled data
  C <- sum(log(choose(n1, y1.1))) + sum(log(choose(n2, y1.2))) -
       sum(log(choose(eba.p$n, eba.p$y1)))  # combinatorial constant

  tests <- rbind(
  # # Discarded, because the remaining deviances do not sum to the
  # #   "overall deviance" associated with this test
  # # mean poisson model vs. saturated poisson model (on y) for both orders
  #  c(df1 <- 1,
  #    df2 <- 2 * I*(I - 1),
  #    l.1 <- sum(dpois(y, mean(y), log=TRUE)),
  #    l.2 <- sum(dpois(y, y, log=TRUE)),
  #    dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # order-effect EBA model vs. saturated binomial model
    c(df1 <- J - 1 + 1,
      df2 <- I*(I - 1),
      l.1 <- x$logL.eba,
      l.2 <- x$logL.sat,
      dev <- 2 * (l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # no-order-effect (= pooled) EBA model vs. order-effect EBA model
    c(df1 <- J - 1,
      df2 <- J - 1 + 1,
      l.1 <- eba.p$logL.eba + C,  # add constant to obtain correct logLik
      l.2 <- x$logL.eba,
      dev <- 2 * (l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # Null model with order effect vs. order-effect EBA model
    c(df1 <- 1,
      df2 <- J -1 + 1,
      l.1 <- sum(dbinom(x$y1, x$n, c(rep(1/(1+order0), npairs),
                 rep(order0/(1+order0), npairs) ), log=TRUE)),
      l.2 <- x$logL.eba,
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1)),

   # mean poisson model vs. saturated poisson model (on n)
    c(df1 <- 1,
      df2 <- I * (I-1),
      l.1 <- sum(dpois(x$n, mean(x$n), log=TRUE)),
      l.2 <- sum(dpois(x$n, x$n, log=TRUE)),
      dev <- 2*(l.2 - l.1), 1 - pchisq(dev, df2 - df1))
  )
  #rownames(tests) <- c("Overall", "EBA.order", "Order", "Effect", "Imbalance")
  rownames(tests) <- c("EBA.order", "Order", "Effect", "Imbalance")
  colnames(tests) <- c("Df1","Df2","logLik1","logLik2","Deviance","Pr(>Chi)")

  aic <- -2*x$logL.eba + 2*(length(coef)-1+1)
  ans <- list(coefficients=coef.table, order.effects=order.table, aic=aic,
    logL.eba=x$logL.eba, logL.sat=x$logL.sat, tests=tests, chi.alt=x$chi.alt)
  class(ans) <- "summary.eba.order"
  return(ans)
}


# Old:
# print.summary.eba.order <- function(x, digits=max(3, getOption("digits")-3),
#   na.print="", symbolic.cor=p>4, signif.stars=getOption("show.signif.stars"),
#   ...){

print.summary.eba.order <- function(x, digits=max(3, getOption("digits") - 3),
  na.print="", signif.stars=getOption("show.signif.stars"), ...){
  cat("\nParameter estimates (H0: parameter = 0):\n")
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
  cat("\nOrder effects (H0: parameter = 1):\n")
  printCoefmat(x$order, digits = digits, signif.stars = signif.stars, ...)
  cat("\nModel tests:\n")
  printCoefmat(x$tests, digits = digits, signif.stars = signif.stars,
    zap.ind = 1:2, tst.ind = 3:5, ...)
  cat("\nAIC: ", format(x$aic,digits=max(4, digits + 1)), "\n")
  cat("Pearson X2:", format(x$chi.alt, digits=digits))
  cat("\n")
  invisible(x)
}
