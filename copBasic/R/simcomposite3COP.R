"simcomposite3COP" <-
function(n=1000, nsim=100, compositor=composite3COP,
         parents=NULL, ploton=FALSE, points=FALSE,
         showpar=FALSE, showresults=FALSE, digits=6,
         ...) {

  if(nsim == Inf) {
    warning("can not have an infinite length (row) matrix")
    return(NA)
  }

  to.get.a.width1 <- parents$para1gen()
  w1 <- length(to.get.a.width1)
  to.get.a.width2 <- parents$para2gen()
  w2 <- length(to.get.a.width2)

  vals <- matrix(nrow=nsim, ncol=(w1 + w2 + 10))
  colnames(vals) <- c("alpha", "beta",
                      "kappa", "gamma",
                      "T2.12", "T2.21",
                      "T3.12", "T3.21",
                      "T4.12", "T4.21",
                      rep("Cop1Thetas",w1),
                      rep("Cop2Thetas",w2))
  i <- 0
  while(i < nsim) {
    Theta1 <- parents$para1gen()
    Theta2 <- parents$para2gen()
    alpha  <- runif(1)
    beta   <- runif(1)
    kappa  <- runif(1)
    gamma  <- runif(1)
    para   <- list(cop1=parents$cop1,
                   cop2=parents$cop2,
                   alpha=alpha,
                   beta=beta,
                   kappa=kappa,
                   gamma=gamma,
                   para1=Theta1,
                   para2=Theta2)

    if(showpar) print(para)

    S <- simCOP(n=n, cop=compositor,
                ploton=ploton, points=points,
                para=para, col=rgb(0,0,0,0.1), pch=16)
    #abline(0,1); abline(1,-1)
    #mtext(paste(c("Theta1 = ",round(Theta1,digits=3), "   ",
    #              "Theta2 = ",round(Theta2,digits=3), "   ",
    #              "Alpha = ", round(alpha,digits=3),  "   ",
    #              "Beta = ",  round(beta,digits=3), "   ",
    #              "Kappa= ",  round(kappa,digits=3), "   ",
    #              "Gamma= ",  round(gamma,digits=3)), sep="",
    #                 collapse=" "))
    #cat(c("Sleeping to continue\n"))

    z <- lmomco::lcomoms2(S, nmom=4, opdiag=TRUE)

    results <- c(alpha, beta, kappa, gamma,
                 z$T2, z$T3, z$T4, Theta1, Theta2)
    if(showresults) cat(c(round(results, digits=digits), "\n"))
    i <- i + 1
    vals[i,] <- results
  }
  return(vals)
}
