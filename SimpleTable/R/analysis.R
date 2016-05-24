##  observed data:
##
##        Y=0    Y=1
## X=0    C00    C01
## X=1    C10    C11
##
## internal cells above are observed frequencies
##
## It is assumed there is an unobserved confounding variable Z that
## corresponds to the following distribution of potential outcomes for Y:
##
##  Y(0)    Y(1)    Z
## -------------------
##   0        0     0
##   0        1     1
##   1        0     2
##   1        1     3
##
##
## Joint distribution for (X,Y):
##
## P(X=0, Y=0) = theta00
## P(X=0, Y=1) = theta01
## P(X=1, Y=0) = theta10
## P(X=1, Y=1) = theta11
##
##
## Conditional distribution for (Z | X, Y):
##
## P(Z=1 | X=0, Y=0) = psi00
## P(Z=3 | X=0, Y=1) = psi01
## P(Z=2 | X=1, Y=0) = psi10
## P(Z=3 | X=1, Y=1) = psi11
## P(Z=0 | X=0, Y=0) = (1-psi00)
## P(Z=2 | X=0, Y=1) = (1-psi01)
## P(Z=0 | X=1, Y=0) = (1-psi10)
## P(Z=1 | X=1, Y=1) = (1-psi11)
##
##
## Priors:
##
## (theta00, theta01, theta10, theta11) ~ Dirichlet(a00, a01, a10, a11)
## psi00 ~ Beta(b00, c00)
## psi01 ~ Beta(b01, c01)
## psi10 ~ Beta(b10, c10)
## psi11 ~ Beta(b11, c11)
##
##
## Interpretatoin of prior for psi:
## Conditional distribution for (Z | X, Y):
##
## Assuming: X=0 placebo, X=1 treatment;  Y=0 fail, Y=1 succeed
##
##  Y(0)    Y(1)    Z
## -------------------
##   0        0     0  Never Succeed
##   0        1     1  Helped
##   1        0     2  Hurt
##   1        1     3  Always Succeed
##
## P(Z=1 | X=0, Y=0) = psi00      b00 Helped
## P(Z=0 | X=0, Y=0) = (1-psi00)  c00 Never Succeed 
##
## P(Z=3 | X=0, Y=1) = psi01      b01 Always Succeed
## P(Z=2 | X=0, Y=1) = (1-psi01)  c01 Hurt
##
## P(Z=2 | X=1, Y=0) = psi10      b10 Hurt
## P(Z=0 | X=1, Y=0) = (1-psi10)  c10 Never Succeed
##
## P(Z=3 | X=1, Y=1) = psi11      b11 Always Succeed
## P(Z=1 | X=1, Y=1) = (1-psi11)  c11 Helped
##
## Monotonicity:
## If very few people are Helped then
##    c00 >> b00 & b11 >> c11
##
## If very few people are Hurt then
##    b01 >> c01 & c10 >> b10
##

"analyze2x2" <- function(C00, C01, C10, C11,
                         a00, a01, a10, a11,
                         b00, b01, b10, b11,
                         c00, c01, c10, c11,
                         nsamp=50000){

  theta.post <- rdirichlet(nsamp, c(C00+a00, C01+a01, C10+a10, C11+a11))
  theta00 <- theta.post[,1]
  theta01 <- theta.post[,2]
  theta10 <- theta.post[,3]
  theta11 <- theta.post[,4]

  theta.prior <- rdirichlet(nsamp, c(a00, a01, a10, a11))
  theta00.prior <- theta.prior[,1]
  theta01.prior <- theta.prior[,2]
  theta10.prior <- theta.prior[,3]
  theta11.prior <- theta.prior[,4]
  
  psi00 <- rbeta(nsamp, b00, c00)
  psi01 <- rbeta(nsamp, b01, c01)
  psi10 <- rbeta(nsamp, b10, c10)
  psi11 <- rbeta(nsamp, b11, c11)


  ## prima facie post intervention distribution for all units
  PostIntDist.pf <- data.frame(PY0.X0=(theta00/(theta00+theta01)),
                               PY1.X0=(theta01/(theta00+theta01)),
                               PY0.X1=(theta10/(theta10+theta11)),
                               PY1.X1=(theta11/(theta10+theta11)))


  ## sensitivity analysis post. int. dist. for all units
  PY1.X1.sens <- theta00*psi00 + theta01*psi01 + theta11
  PY0.X1.sens <-  1 - PY1.X1.sens
  PY1.X0.sens <- theta10*psi10 + theta11*psi11 + theta01
  PY0.X0.sens <- 1 - PY1.X0.sens

  PostIntDist.sens <- data.frame(PY0.X0=PY0.X0.sens,
                                 PY1.X0=PY1.X0.sens,
                                 PY0.X1=PY0.X1.sens,
                                 PY1.X1=PY1.X1.sens)


  ## sensitivity analysis post. int. dist. for treated units
  PY1.X1.sens <- theta11 / (theta10 + theta11)
  PY0.X1.sens <- 1 - PY1.X1.sens
  PY1.X0.sens <- (psi10*theta10 + psi11*theta11) / (theta10 + theta11)
  PY0.X0.sens <- 1 - PY1.X1.sens
  
  PostIntDist.T.sens <- data.frame(PY0.X0=PY0.X0.sens,
                                   PY1.X0=PY1.X0.sens,
                                   PY0.X1=PY0.X1.sens,
                                   PY1.X1=PY1.X1.sens)
  
  
  ## sensitivity analysis post. int. dist. for control units
  PY1.X1.sens <- (psi00*theta00 + psi01*theta01) / (theta00 + theta01)
  PY0.X1.sens <- 1 - PY1.X1.sens
  PY1.X0.sens <- theta01 / (theta00 + theta01)
  PY0.X0.sens <- 1 - PY1.X1.sens
  
  PostIntDist.C.sens <- data.frame(PY0.X0=PY0.X0.sens,
                                   PY1.X0=PY1.X0.sens,
                                   PY0.X1=PY0.X1.sens,
                                   PY1.X1=PY1.X1.sens)
  
  


  ## prior post. int. dist. for all units
  PY1.X1.prior <- theta00.prior*psi00 + theta01.prior*psi01 + theta11.prior
  PY0.X1.prior <-  1 - PY1.X1.prior
  PY1.X0.prior <- theta10.prior*psi10 + theta11.prior*psi11 + theta01.prior
  PY0.X0.prior <- 1 - PY1.X0.prior

  PostIntDist.prior <- data.frame(PY0.X0=PY0.X0.prior,
                                  PY1.X0=PY1.X0.prior,
                                  PY0.X1=PY0.X1.prior,
                                  PY1.X1=PY1.X1.prior)
  

  ## prior post. int. dist. for treated units
  PY1.X1.prior <- theta11.prior / (theta10.prior + theta11.prior)
  PY0.X1.prior <- 1 - PY1.X1.prior
  PY1.X0.prior <- (psi10*theta10.prior + psi11*theta11.prior) /
    (theta10.prior + theta11.prior)
  PY0.X0.prior <- 1 - PY1.X1.prior
  
  PostIntDist.T.prior <- data.frame(PY0.X0=PY0.X0.prior,
                                    PY1.X0=PY1.X0.prior,
                                    PY0.X1=PY0.X1.prior,
                                    PY1.X1=PY1.X1.prior)
  
  
  ## prior post. int. dist. for control units
  PY1.X1.prior <- (psi00*theta00.prior + psi01*theta01.prior) /
    (theta00.prior + theta01.prior)
  PY0.X1.prior <- 1 - PY1.X1.prior
  PY1.X0.prior <- theta01.prior / (theta00.prior + theta01.prior)
  PY0.X0.prior <- 1 - PY1.X1.prior
  
  PostIntDist.C.prior <- data.frame(PY0.X0=PY0.X0.prior,
                                    PY1.X0=PY1.X0.prior,
                                    PY0.X1=PY0.X1.prior,
                                    PY1.X1=PY1.X1.prior)
  

  
  
  ## calculate large sample nonparametric bounds
  n <- sum(C00, C01, C10, C11)
  ATE.min <- -(C01/n + C10/n)
  ATE.max <- C00/n + C11/n
  ATT.min <- - C10 / (C10 + C11)
  ATT.max <- C11 / (C10 + C11)
  ATC.min <- -C01 / (C00 + C01)
  ATC.max <- C00 / (C00 + C01)
  
  RR.min <- C11 / (C01 + C10 + C11)
  RR.max <- (C11 + C00 + C01) / C01
  RRT.min <- C11 / (C10 + C11)
  RRT.max <- Inf
  RRC.min <- 0
  RRC.max <- 1 + C00 / C01
  

  

  ## calculate causal effects
  ## prima facie effects
  ATE.pf <- PostIntDist.pf$PY1.X1 - PostIntDist.pf$PY1.X0
  ATT.pf <- ATE.pf
  ATC.pf <- ATE.pf
  RR.pf <- PostIntDist.pf$PY1.X1 / PostIntDist.pf$PY1.X0
  RRT.pf <- RR.pf
  RRC.pf <- RR.pf


  ## sensitivity analysis effects
  ATE.sens <- PostIntDist.sens$PY1.X1 - PostIntDist.sens$PY1.X0
  ATT.sens <- PostIntDist.T.sens$PY1.X1 - PostIntDist.T.sens$PY1.X0
  ATC.sens <- PostIntDist.C.sens$PY1.X1 - PostIntDist.C.sens$PY1.X0
  RR.sens <- PostIntDist.sens$PY1.X1 / PostIntDist.sens$PY1.X0
  RRT.sens <- PostIntDist.T.sens$PY1.X1 / PostIntDist.T.sens$PY1.X0
  RRC.sens <- PostIntDist.C.sens$PY1.X1 / PostIntDist.C.sens$PY1.X0
  

  ## prior effects
  ATE.prior <- PostIntDist.prior$PY1.X1 - PostIntDist.prior$PY1.X0
  ATT.prior <- PostIntDist.T.prior$PY1.X1 - PostIntDist.T.prior$PY1.X0
  ATC.prior <- PostIntDist.C.prior$PY1.X1 - PostIntDist.C.prior$PY1.X0
  RR.prior <- PostIntDist.prior$PY1.X1 / PostIntDist.prior$PY1.X0
  RRT.prior <- PostIntDist.T.prior$PY1.X1 / PostIntDist.T.prior$PY1.X0
  RRC.prior <- PostIntDist.C.prior$PY1.X1 / PostIntDist.C.prior$PY1.X0
  
  
  return(structure(list(C00=C00, C01=C01, C10=C10, C11=C11,
                        a00=a00, a01=a01, a10=a10, a11=a11,
                        b00=b00, b01=b01, b10=b10, b11=b11,
                        c00=c00, c01=c01, c10=c10, c11=c11,
                        theta00=theta00, theta01=theta01,
                        theta10=theta10, theta11=theta11,
                        theta00.prior=theta00.prior,
                        theta01.prior=theta01.prior,
                        theta10.prior=theta10.prior,
                        theta11.prior=theta11.prior,
                        psi00=psi00, psi01=psi01, psi10=psi10, psi11=psi11,
                        PostIntDist.pf=PostIntDist.pf,
                        PostIntDist.sens=PostIntDist.sens,
                        PostIntDist.T.sens=PostIntDist.T.sens,
                        PostIntDist.C.sens=PostIntDist.C.sens,
                        PostIntDist.prior=PostIntDist.prior,
                        PostIntDist.T.prior=PostIntDist.T.prior,
                        PostIntDist.C.prior=PostIntDist.C.prior,
                        ATE.min=ATE.min,
                        ATE.max=ATE.max,
                        ATT.min=ATT.min,
                        ATT.max=ATT.max,
                        ATC.min=ATC.min,
                        ATC.max=ATC.max,
                        RR.min=RR.min,
                        RR.max=RR.max,
                        RRT.min=RRT.min,
                        RRT.max=RRT.max,
                        RRC.min=RRC.min,
                        RRC.max=RRC.max,
                        logRR.min=log(RR.min),
                        logRR.max=log(RR.max),
                        logRRT.min=log(RRT.min),
                        logRRT.max=log(RRT.max),
                        logRRC.min=log(RRC.min),
                        logRRC.max=log(RRC.max),                        
                        ATE.pf=ATE.pf,
                        ATT.pf=ATT.pf,
                        ATC.pf=ATC.pf,
                        RR.pf=RR.pf,
                        RRT.pf=RRT.pf,
                        RRC.pf=RRC.pf,
                        logRR.pf=log(RR.pf),
                        logRRT.pf=log(RRT.pf),
                        logRRC.pf=log(RRC.pf),
                        ATE.sens=ATE.sens,
                        ATT.sens=ATT.sens,
                        ATC.sens=ATC.sens,
                        RR.sens=RR.sens,
                        RRT.sens=RRT.sens,
                        RRC.sens=RRC.sens,
                        logRR.sens=log(RR.sens),
                        logRRT.sens=log(RRT.sens),
                        logRRC.sens=log(RRC.sens),
                        ATE.prior=ATE.prior,
                        ATT.prior=ATT.prior,
                        ATC.prior=ATC.prior,
                        RR.prior=RR.prior,
                        RRT.prior=RRT.prior,
                        RRC.prior=RRC.prior,
                        logRR.prior=log(RR.prior),
                        logRRT.prior=log(RRT.prior),
                        logRRC.prior=log(RRC.prior)
                        ),
                   class="SimpleTable"))
  
  
} ## end analyze2x2











## function to deal with situation where there is a single
## measured categorical confounder (W) in addition to Z
##
## works by combining multiple SimpleTable objects created
## from (X,Y,Z|W=w) across all levels of W
analyze2x2xK <- function(SimpleTableList, Wpriorvector){

  K <- length(SimpleTableList)
  if (K != length(Wpriorvector)){
    stop("SimpleTableList and Wpriorvector not of same length.\n")
  }

  Wcountvec <- rep(NA, K)
  for (w in 1:K){
    Wcountvec[w] <- SimpleTableList[[w]]$C00 + SimpleTableList[[w]]$C01 +
      SimpleTableList[[w]]$C10 + SimpleTableList[[w]]$C11
  }

  ntotal <- sum(Wcountvec)

  ## use the smallest number of samples among the stratum specific
  ## results as the number of samples of psi (the marginal probability
  ## of W
  nsamp <- length(SimpleTableList[[1]]$theta00)
  for (j in 2:K){
    nsamp <- min(length(SimpleTableList[[j]]$theta00), nsamp)
  }

  ## W is asssumed multinomial with a Dirichlet prior
  ## the following samples from the posterior
  phi.mat <- rdirichlet(nsamp, Wcountvec+Wpriorvector)
  phi.mat.prior <- rdirichlet(nsamp, Wpriorvector)

  PY0.X0.pf <- PY0.X0.sens <- PY0.X0.T.sens <- PY0.X0.C.sens <-
    PY0.X0.prior <- PY0.X0.T.prior <- PY0.X0.C.prior <- 0
  PY1.X0.pf <- PY1.X0.sens <- PY1.X0.T.sens <- PY1.X0.C.sens <-
    PY1.X0.prior <- PY1.X0.T.prior <- PY1.X0.C.prior <- 0
  PY0.X1.pf <- PY0.X1.sens <- PY0.X1.T.sens <- PY0.X1.C.sens <-
    PY0.X1.prior <- PY0.X1.T.prior <- PY0.X1.C.prior <- 0
  PY1.X1.pf <- PY1.X1.sens <- PY1.X1.T.sens <- PY1.X1.C.sens <-
    PY1.X1.prior <- PY1.X1.T.prior <- PY1.X1.C.prior <- 0
  C00 <- C01 <- C10 <- C11 <- 0

  ## average over all K of the analyses conditional on levels of W
  for (w in 1:K){

    PY1.X0.pf <- PY1.X0.pf + (SimpleTableList[[w]]$theta01[1:nsamp] /
                              (SimpleTableList[[w]]$theta00[1:nsamp] +
                               SimpleTableList[[w]]$theta01[1:nsamp])) *
                                       phi.mat[,w]

    PY1.X1.pf <- PY1.X1.pf + (SimpleTableList[[w]]$theta11[1:nsamp] /
                              (SimpleTableList[[w]]$theta10[1:nsamp] +
                               SimpleTableList[[w]]$theta11[1:nsamp])) *
                                 phi.mat[,w]



    PY1.X1.sens <- PY1.X1.sens + (SimpleTableList[[w]]$theta00[1:nsamp] *
                                  SimpleTableList[[w]]$psi00[1:nsamp] +
                                  SimpleTableList[[w]]$theta01[1:nsamp] *
                                  SimpleTableList[[w]]$psi01[1:nsamp] +
                                  SimpleTableList[[w]]$theta11[1:nsamp]) *
                                    phi.mat[,w]
    
    PY1.X0.sens <- PY1.X0.sens + (SimpleTableList[[w]]$theta10[1:nsamp] *
                                  SimpleTableList[[w]]$psi10[1:nsamp] +
                                  SimpleTableList[[w]]$theta11[1:nsamp] *
                                  SimpleTableList[[w]]$psi11[1:nsamp] +
                                  SimpleTableList[[w]]$theta01[1:nsamp]) *
                                    phi.mat[,w]

    PY1.X1.T.sens <- PY1.X1.T.sens + (SimpleTableList[[w]]$theta11[1:nsamp] /
                                      (SimpleTableList[[w]]$theta10[1:nsamp] +
                                       SimpleTableList[[w]]$theta11[1:nsamp])) *
                                         phi.mat[,w]
    
    PY1.X0.T.sens <- PY1.X0.T.sens + ((SimpleTableList[[w]]$psi10 *
                                       SimpleTableList[[w]]$theta10 +
                                       SimpleTableList[[w]]$psi11 *
                                       SimpleTableList[[w]]$theta11) /
                                      (SimpleTableList[[w]]$theta10 +
                                       SimpleTableList[[w]]$theta11)) *
                                         phi.mat[,w]
    
    PY1.X1.C.sens <- PY1.X1.C.sens + ((SimpleTableList[[w]]$psi00[1:nsamp] *
                                       SimpleTableList[[w]]$theta00[1:nsamp] +
                                       SimpleTableList[[w]]$psi01[1:nsamp] *
                                       SimpleTableList[[w]]$theta01[1:nsamp]) /
                                      (SimpleTableList[[w]]$theta00[1:nsamp] +
                                       SimpleTableList[[w]]$theta01[1:nsamp])) *
                                         phi.mat[,w]
    
    PY1.X0.C.sens <- PY1.X0.C.sens + (SimpleTableList[[w]]$theta01[1:nsamp] /
                                      (SimpleTableList[[w]]$theta00[1:nsamp] +
                                       SimpleTableList[[w]]$theta01[1:nsamp])) *
                                         phi.mat[,w]

    


    PY1.X1.prior <- PY1.X1.prior+(SimpleTableList[[w]]$theta00.prior[1:nsamp] *
                                  SimpleTableList[[w]]$psi00[1:nsamp] +
                                  SimpleTableList[[w]]$theta01.prior[1:nsamp] *
                                  SimpleTableList[[w]]$psi01[1:nsamp] +
                                  SimpleTableList[[w]]$theta11.prior[1:nsamp]) *
                                    phi.mat.prior[,w]
    
    PY1.X0.prior <- PY1.X0.prior+(SimpleTableList[[w]]$theta10.prior[1:nsamp] *
                                  SimpleTableList[[w]]$psi10[1:nsamp] +
                                  SimpleTableList[[w]]$theta11.prior[1:nsamp] *
                                  SimpleTableList[[w]]$psi11[1:nsamp] +
                                  SimpleTableList[[w]]$theta01.prior[1:nsamp]) *
                                    phi.mat.prior[,w]

    PY1.X1.T.prior<- PY1.X1.T.prior +
      (SimpleTableList[[w]]$theta11.prior[1:nsamp] /
       (SimpleTableList[[w]]$theta10.prior[1:nsamp] +
        SimpleTableList[[w]]$theta11.prior[1:nsamp])) *
          phi.mat.prior[,w]
    
    PY1.X0.T.prior <- PY1.X0.T.prior +
      ((SimpleTableList[[w]]$psi10 *
        SimpleTableList[[w]]$theta10.prior +
        SimpleTableList[[w]]$psi11 *
        SimpleTableList[[w]]$theta11.prior) /
       (SimpleTableList[[w]]$theta10.prior +
        SimpleTableList[[w]]$theta11.prior)) *
          phi.mat.prior[,w]
    
    PY1.X1.C.prior <- PY1.X1.C.prior +
      ((SimpleTableList[[w]]$psi00[1:nsamp] *
        SimpleTableList[[w]]$theta00.prior[1:nsamp] +
        SimpleTableList[[w]]$psi01[1:nsamp] *
        SimpleTableList[[w]]$theta01.prior[1:nsamp]) /
       (SimpleTableList[[w]]$theta00.prior[1:nsamp] +
        SimpleTableList[[w]]$theta01.prior[1:nsamp])) *
          phi.mat.prior[,w]
    
    PY1.X0.C.prior <- PY1.X0.C.prior +
      (SimpleTableList[[w]]$theta01.prior[1:nsamp] /
       (SimpleTableList[[w]]$theta00.prior[1:nsamp] +
        SimpleTableList[[w]]$theta01.prior[1:nsamp])) *
          phi.mat.prior[,w]
    

    
    
    C00 <- C00 + SimpleTableList[[w]]$C00
    C01 <- C01 + SimpleTableList[[w]]$C01
    C10 <- C10 + SimpleTableList[[w]]$C10
    C11 <- C11 + SimpleTableList[[w]]$C11

  } ## end w loop


  PY0.X0.pf <- 1 - PY1.X0.pf
  PY0.X1.pf <- 1 - PY1.X1.pf
  
  PY0.X0.sens <- 1 - PY1.X0.sens
  PY0.X1.sens <- 1 - PY1.X1.sens
  PY0.X0.T.sens <- 1 - PY1.X0.T.sens
  PY0.X1.T.sens <- 1 - PY1.X1.T.sens
  PY0.X0.C.sens <- 1 - PY1.X0.C.sens
  PY0.X1.C.sens <- 1 - PY1.X1.C.sens
  
  PY0.X0.prior <- 1 - PY1.X0.prior
  PY0.X1.prior <- 1 - PY1.X1.prior
  PY0.X0.T.prior <- 1 - PY1.X0.T.prior
  PY0.X1.T.prior <- 1 - PY1.X1.T.prior
  PY0.X0.C.prior <- 1 - PY1.X0.C.prior
  PY0.X1.C.prior <- 1 - PY1.X1.C.prior



  ## prima facie post intervention distribution for all units
  PostIntDist.pf <- data.frame(PY0.X0=PY0.X0.pf,
                               PY1.X0=PY1.X0.pf,
                               PY0.X1=PY0.X1.pf,
                               PY1.X1=PY1.X1.pf)


  ## sensitivity analysis post. int. dist. for all units
  PostIntDist.sens <- data.frame(PY0.X0=PY0.X0.sens,
                                 PY1.X0=PY1.X0.sens,
                                 PY0.X1=PY0.X1.sens,
                                 PY1.X1=PY1.X1.sens)


  ## sensitivity analysis post. int. dist. for treated units
  PostIntDist.T.sens <- data.frame(PY0.X0=PY0.X0.T.sens,
                                   PY1.X0=PY1.X0.T.sens,
                                   PY0.X1=PY0.X1.T.sens,
                                   PY1.X1=PY1.X1.T.sens)
  
  
  ## sensitivity analysis post. int. dist. for control units
  PostIntDist.C.sens <- data.frame(PY0.X0=PY0.X0.C.sens,
                                   PY1.X0=PY1.X0.C.sens,
                                   PY0.X1=PY0.X1.C.sens,
                                   PY1.X1=PY1.X1.C.sens)

  

  ## prior post. int. dist. for all units
  PostIntDist.prior <- data.frame(PY0.X0=PY0.X0.prior,
                                  PY1.X0=PY1.X0.prior,
                                  PY0.X1=PY0.X1.prior,
                                  PY1.X1=PY1.X1.prior)
  

  ## prior post. int. dist. for treated units
  PostIntDist.T.prior <- data.frame(PY0.X0=PY0.X0.T.prior,
                                    PY1.X0=PY1.X0.T.prior,
                                    PY0.X1=PY0.X1.T.prior,
                                    PY1.X1=PY1.X1.T.prior)
  
  
  ## prior post. int. dist. for control units
  PostIntDist.C.prior <- data.frame(PY0.X0=PY0.X0.C.prior,
                                    PY1.X0=PY1.X0.C.prior,
                                    PY0.X1=PY0.X1.C.prior,
                                    PY1.X1=PY1.X1.C.prior)
  

  
  ## calculate large sample nonparametric bounds
  n <- sum(C00, C01, C10, C11)
  ATE.min <- -(C01/n + C10/n)
  ATE.max <- C00/n + C11/n
  ATT.min <- - C10 / (C10 + C11)
  ATT.max <- C11 / (C10 + C11)
  ATC.min <- -C01 / (C00 + C01)
  ATC.max <- C00 / (C00 + C01)
  
  RR.min <- C11 / (C01 + C10 + C11)
  RR.max <- (C11 + C00 + C01) / C01
  RRT.min <- C11 / (C10 + C11)
  RRT.max <- Inf
  RRC.min <- 0
  RRC.max <- 1 + C00 / C01




  ## calculate causal effects
  ## prima facie effects
  ATE.pf <- PostIntDist.pf$PY1.X1 - PostIntDist.pf$PY1.X0
  ATT.pf <- ATE.pf
  ATC.pf <- ATE.pf
  RR.pf <- PostIntDist.pf$PY1.X1 / PostIntDist.pf$PY1.X0
  RRT.pf <- RR.pf
  RRC.pf <- RR.pf


  ## sensitivity analysis effects
  ATE.sens <- PostIntDist.sens$PY1.X1 - PostIntDist.sens$PY1.X0
  ATT.sens <- PostIntDist.T.sens$PY1.X1 - PostIntDist.T.sens$PY1.X0
  ATC.sens <- PostIntDist.C.sens$PY1.X1 - PostIntDist.C.sens$PY1.X0
  RR.sens <- PostIntDist.sens$PY1.X1 / PostIntDist.sens$PY1.X0
  RRT.sens <- PostIntDist.T.sens$PY1.X1 / PostIntDist.T.sens$PY1.X0
  RRC.sens <- PostIntDist.C.sens$PY1.X1 / PostIntDist.C.sens$PY1.X0
  

  ## prior effects
  ATE.prior <- PostIntDist.prior$PY1.X1 - PostIntDist.prior$PY1.X0
  ATT.prior <- PostIntDist.T.prior$PY1.X1 - PostIntDist.T.prior$PY1.X0
  ATC.prior <- PostIntDist.C.prior$PY1.X1 - PostIntDist.C.prior$PY1.X0
  RR.prior <- PostIntDist.prior$PY1.X1 / PostIntDist.prior$PY1.X0
  RRT.prior <- PostIntDist.T.prior$PY1.X1 / PostIntDist.T.prior$PY1.X0
  RRC.prior <- PostIntDist.C.prior$PY1.X1 / PostIntDist.C.prior$PY1.X0




  return(structure(list(PostIntDist.pf=PostIntDist.pf,
                        PostIntDist.sens=PostIntDist.sens,
                        PostIntDist.T.sens=PostIntDist.T.sens,
                        PostIntDist.C.sens=PostIntDist.C.sens,
                        PostIntDist.prior=PostIntDist.prior,
                        PostIntDist.T.prior=PostIntDist.T.prior,
                        PostIntDist.C.prior=PostIntDist.C.prior,
                        ATE.min=ATE.min,
                        ATE.max=ATE.max,
                        ATT.min=ATT.min,
                        ATT.max=ATT.max,
                        ATC.min=ATC.min,
                        ATC.max=ATC.max,
                        RR.min=RR.min,
                        RR.max=RR.max,
                        RRT.min=RRT.min,
                        RRT.max=RRT.max,
                        RRC.min=RRC.min,
                        RRC.max=RRC.max,
                        logRR.min=log(RR.min),
                        logRR.max=log(RR.max),
                        logRRT.min=log(RRT.min),
                        logRRT.max=log(RRT.max),
                        logRRC.min=log(RRC.min),
                        logRRC.max=log(RRC.max),                        
                        ATE.pf=ATE.pf,
                        ATT.pf=ATT.pf,
                        ATC.pf=ATC.pf,
                        RR.pf=RR.pf,
                        RRT.pf=RRT.pf,
                        RRC.pf=RRC.pf,
                        logRR.pf=log(RR.pf),
                        logRRT.pf=log(RRT.pf),
                        logRRC.pf=log(RRC.pf),
                        ATE.sens=ATE.sens,
                        ATT.sens=ATT.sens,
                        ATC.sens=ATC.sens,
                        RR.sens=RR.sens,
                        RRT.sens=RRT.sens,
                        RRC.sens=RRC.sens,
                        logRR.sens=log(RR.sens),
                        logRRT.sens=log(RRT.sens),
                        logRRC.sens=log(RRC.sens),
                        ATE.prior=ATE.prior,
                        ATT.prior=ATT.prior,
                        ATC.prior=ATC.prior,
                        RR.prior=RR.prior,
                        RRT.prior=RRT.prior,
                        RRC.prior=RRC.prior,
                        logRR.prior=log(RR.prior),
                        logRRT.prior=log(RRT.prior),
                        logRRC.prior=log(RRC.prior)
                        ),
                   class="SimpleTable"))

  
  
} ## end analyze2x2xK










"is.SimpleTable" <- function(S){
  return(class(S) == "SimpleTable")
}







