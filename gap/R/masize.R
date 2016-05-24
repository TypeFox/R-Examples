antilogit <- function(x) 1/(1+exp(-x))

getPTE <- function(b1, b2, rho, sdx1=1, sdx2=1) b2*sdx2*rho/((b1+b2*sdx2*rho/sdx1)*sdx1)

getb1star <- function(b1, b2, rho, sdx1=1, sdx2=1) b1+b2*sdx2*rho/sdx1

getb0 <- function(p, X, b1, b2)
{
      b0 <- 1
      p.test <- 0.9999
      while(p.test>p) {
            b0 <- b0-0.1
            p.test <- mean(antilogit(X%*%as.matrix(c(b0, b1, b2))))
        }
      while(p.test<p) {
            b0 <- b0+0.01
            p.test <- mean(antilogit(X%*%as.matrix(c(b0, b1, b2))))
        }
      while(p.test>p) {
            b0 <- b0-0.001
            p.test <- mean(antilogit(X%*%as.matrix(c(b0, b1, b2))))
        }
      while(p.test<p) {
            b0 <- b0+0.0001
            p.test <- mean(antilogit(X%*%as.matrix(c(b0, b1, b2))))
        }
      b0
}

masize <- function(model,opts, alpha=0.025, gamma=0.2)
{
   for(p in c("survival")) {
      if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
         if (!require(p, quietly = TRUE, character.only=TRUE))
         warning(paste("masize needs package `", p, "' to be fully functional; please install", sep=""))
      }
   }
   linear <- paste("linear",1:4,sep="")
   logistic <- paste("logistic",1:5,sep="")
   poisson <- paste("poisson",1:9,sep="")
   cox <- paste("cox",1:9,sep="")
   model.int <- charmatch(model, c(linear,logistic,poisson,cox))
   if (is.na(model.int)) stop("Invalid model type")
   if (model.int == 0) stop("Ambiguous model type")
   if (!is.null(opts$seed)) set.seed(opts$seed)
   ## linear models
   if (model.int == 1) # (2) linear
   {
       b2 <- opts$b2       # regression coefficient for mediator in full model
       rho <- opts$rho     # correlation of primary predictor and mediator
       sdx2 <- opts$sdx2   # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
       sdy <- opts$sdy     # SD of outcome
       desc <- "linear"
       n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2*(1-rho^2))
   }
   if (model.int == 2) # (4) lineara
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        PTE <- opts$PTE       # proportion of effect of primary predictor explained by mediator: (b1star - b1)/b1star
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f
        sdy <- opts$sdy       # SD of outcome
        desc <- "lineara"
        n <- (qnorm(alpha)+qnorm(gamma))^2*rho^2*sdy^2/((b1star*sdx1*PTE)^2*(1-rho^2))
   }
   if (model.int == 3) # (5) linearb
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        b2 <- opts$b2         # regression coefficient for mediator in full model
        PTE <- opts$PTE       # proportion of effect of primary predictor explained by mediator: (b1star - b1)/b1star
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f1
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        sdy <- opts$sdy       # SD of outcome
        desc <- "linearb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2-(b1star*sdx1*PTE)^2)
   }
   if (model.int == 4) # (6) linearc
   {
        b1star <- opts$b1star # regression coefficient for primary predictor in reduced model
        b1 <- opts$b1         # regression coefficient for primary predictor in full model
        b2 <- opts$b2         # regression coefficient for mediator in full model
        sdx1 <- opts$sdx1     # SD of primary predictor; use sqrt(f1*(1-f1)) for binary predictor with prevalence f1
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        sdy <- opts$sdy       # SD of outcome
        desc <- "linearc"
        n <- (qnorm(alpha)+qnorm(gamma))^2*sdy^2/((b2*sdx2)^2-((b1star-b1)*sdx1)^2)
   }
   ## logistic model
   if (model.int == 5) # (7) logistic.approx
   {
        p <- opts$p           # marginal prevalence of outcome
        b2 <- opts$b2         # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2     # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "logistic.approx"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*p*(1-p))
   }
   if (model.int == 6) # (8) logistic.ccs
   {
        p <- opts$p           # marginal prevalence of outcome
        b1 <- opts$b1         # regression coefficient (log odds-ratio) for primary predictor
        b2 <- opts$b2         # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho       # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1     # SD of primary predictor
        sdx2 <- opts$sdx2     # SD of mediator
        ns <- opts$ns         # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- cbind(rep(1, ns), matrix(rnorm(2*ns), ns, 2) %*% chol(S))
        b0 <- getb0(p, X, b1, b2)
        pi <- antilogit(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.ccs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 7) # (8) logistic.bcs
   {
        p <- opts$p                      # marginal prevalence of outcome
        b1 <- opts$b1                    # regression coefficient (log odds-ratio) for primary predictor
        f1 <- opts$f1                    # prevalence of binary primary predictor
        b2 <- opts$b2                    # regression coefficient (log odds-ratio) for mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2                # SD of mediator
        ns <- opts$ns                    # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1)) # conditional mean of mediator when X1=0
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)  # conditional mean of mediator when X1=1
        sdx2.1 <- sdx2*sqrt(1-rho^2)     # SD of mediator within levels of X1
        n0 <- round(ns*(1-f1))           # number of simulated observations with X1=0
        n1 <- ns-n0                      # number of simulated observations with X1=1
        X <- rbind(cbind(rep(1, n0), rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        b0 <- getb0(p, X, b1, b2)
        pi <- antilogit(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.bcs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   if (model.int == 8) # (8) logistic.cbs
   {
        p <- opts$p                     # marginal prevalence of outcome
        b1 <- opts$b1                   # regression coefficient (log odds-ratio) for continuous primary predictor
        b2 <- opts$b2                   # regression coefficient (log odds-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of binary mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1               # SD of continuous primary predictor
        ns <- opts$ns                   # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n0 <- round(ns*(1-f2))
        n1 <- ns-n0
        X <- rbind(cbind(rep(1, n0), rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        b0 <- getb0(p, X, b1, b2)
        pi <- antilogit(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.cbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   if (model.int == 9) # (8) logistic.bbs
   {
        p <- opts$p                     # marginal prevalence of outcome
        b1 <- opts$b1                   # regression coefficient (log odds-ratio) for binary primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log odds-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        ns <- opts$ns                   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2 # fraction with X1=1 and X2=1
        f10 <- f1-f11                   # fraction with X1=1 and X2=0
        f01 <- f2-f11                   # fraction with X1=0 and X2=1
        n11 <- round(ns*f11)
        n10 <- round(ns*f10)
        n01 <- round(ns*f01)
        n00 <- max(1, ns-n10-n01-n11)
        X <- rbind(cbind(rep(1, n00), rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(1, n10), rep(0, n10)),
                   cbind(rep(1, n01), rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11), rep(1, n11)))
        b0 <- getb0(p, X, b1, b2)
        pi <- antilogit(X%*%as.matrix(c(b0, b1, b2)))
        XVX <- crossprod(X*as.vector(sqrt(pi*(1-pi))))
        desc <- "logistic.bbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/b2^2
   }
   ## Poisson model
   if (model.int == 10) # (9) poisson.approx
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.approx"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 11) # (9) poisson.cc
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.cc"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 12) # (9) poisson.bc
   {
        m <- opts$m                     # marginal mean of outcome
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator, use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "poisson.bc"
        n <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2)*m)
   }
   if (model.int == 13) # (10) poisson.cb
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for binary mediator
        f2 <- opts$f2                    # prevalence of binary mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx1 <- opts$sdx1                # SD of continuous primary predictor
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2)) # mean of X1 when X2 = 0
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)  # mean of X1 when X2 = 1
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        B <- f2*exp(b1*mu1+b2)
        D <- (1-f2)*exp(b1*mu0)
        desc <- "poisson.cb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*(((B+D)*sdx1.2)^2+B*D*(mu1-mu0)^2)/(b2^2*B*D*sdx1.2^2*m)
   }
   if (model.int == 14) # (11) poisson.bb
   {
        m <- opts$m                     # marginal mean of outcome
        b1 <- opts$b1                   # regression coefficient (log rate-ratio) for binary primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for binary mediator
        f2 <- opts$f2                   # prevalence of binary mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        f11 <- f1*f2+rho*sqrt(f1*(1-f1)*f2*(1-f2))
        f10 <- f1-f11
        f01 <- f2-f11
        f00 <- 1-f01-f10-f11
        B <- f00
        C <- f10*exp(b1)
        D <- f01*exp(b2)
        E <- f11*exp(b1+b2)
        desc <- "poisson.bb"
        n <- (qnorm(alpha)+qnorm(gamma))^2*(B+C+D+E)*(B+D)*(C+E)/(b2^2*(B*C*D+B*C*E+B*D*E+C*D*E)*m)
   }
   if (model.int == 15) # (12) poisson.ccs
   {
        m <- opts$m     # marginal mean of outcome
        b1 <- opts$b1   # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2   # regression coefficient (log rate-ratio) for continuous mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        sdx1 <- opts$adx1 # SD of primary predictor
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.ccs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 16) # (12) poisson.bcs
   {
        m <- opts$m                     # marginal mean of outcome
        b1 <- opts$b1                   # regression coefficient (log rate-ratio) for continuous primary predictor
        f1 <- opts$f1                   # prevalence of primary predictor
        b2 <- opts$b2                   # regression coefficient (log rate-ratio) for continuous mediator
        rho <- opts$rho                 # correlation of primary predictor and mediator
        sdx2 <- opts$sdx2               # SD of mediator
        ns <- opts$ns                   # number of observations in simulated dataset
        n1 <- round(ns*f1)
        n0 <- ns-n1
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1)) # mean of X2 when X1 = 0
        mu1 <- rho*sdx2*sqrt((1-f1)/f1) # mean of X2 when X1 = 1
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.bcs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 17) # (12) poisson.cbs
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for continuous primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for continuous mediator
        f2 <- opts$f2                    # prevalence of mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        sdx1 <- opts$sdc1                # SD of primary predictor
        ns <- opts$ns                    # number of observations in simulated dataset
        n1 <- round(ns*f2)
        n0 <- ns-n1
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2)) # mean of X1 when X2 = 0
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)  # mean of X1 when X2 = 1
        sdx1.2 = sdx1*sqrt(1-rho^2)
        X <-  rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                    cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.cbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 18) # (12) poisson.bbs
   {
        m <- opts$m                      # marginal mean of outcome
        b1 <- opts$b1                    # regression coefficient (log rate-ratio) for primary predictor
        f1 <- opts$f1                    # prevalence of primary predictor
        b2 <- opts$b2                    # regression coefficient (log rate-ratio) for mediator
        f2 <- opts$f2                    # prevalence of mediator
        rho <- opts$rho                  # correlation of primary predictor and mediator
        ns <- opts$ns                    # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        f10 <- f1-f11
        f01 <- f2-f11
        f00 <- 1-f01-f10-f11
        n11 <- round(f11*ns)
        n01 <- round(f01*ns)
        n10 <- round(f10*ns)
        n00 <- ns-n01-n10-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        rr <- exp(X%*%matrix(c(b1, b2)))
        XVX <- crossprod(cbind(rep(1, ns), X)*as.vector(sqrt(m*rr/mean(rr))))
        desc <- "poisson.bbs"
        n <- (qnorm(alpha)+qnorm(gamma))^2*solve(XVX/ns)[3,3]/(b2^2)
   }
   if (model.int == 19) # (13) cox.approx
   {
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx2 <- opts$sdx2 # SD of mediator; use sqrt(f2*(1-f2)) for binary mediator with prevalence f2
        desc <- "cox.approx"
        d <- (qnorm(alpha)+qnorm(gamma))^2/((b2*sdx2)^2*(1-rho^2))
        n <- d/f
   }
   if (model.int == 20) # (14) cox.ccs
   {
        b1 <- opts$b1      # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2      # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho    # correlation of primary predictor and mediator
        f <- opts$f        # fraction of follow-up times that are uncensored
        sdx1 <- opts$sdx1  # SD of primary predictor
        sdx2 <- opts$sdx2  # SD of mediator
        ns <- opts$ns      # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.ccs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 21) # (14) cox.bcs
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1     # prevalence of primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1))
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        n1 <- round(f1*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bcs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 22) # (14) cox.cbs
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2     # prevalence of mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        sdx1 <- opts$sdx1 # SD of primary predictor
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n1 <- round(f2*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.cbs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 23) # (14) cox.bbs
   {
        b1 <- opts$b1   # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1   # prevalence of primary predictor
        b2 <- opts$b2   # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2   # prevalence of mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        f <- opts$f     # fraction of follow-up times that are uncensored
        ns <- opts$ns   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        n11 <- round(ns*f11)
        n10 <- round(ns*(f1-f11))
        n01 <- round(ns*(f2-f11))
        n00 <- ns-n10-n01-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        time <- rexp(ns, rate = exp(X%*%as.matrix(c(b1, b2))))
        event <- ifelse(rank(time)<=round(ns*f), 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bbs"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 24) # (14) cox.ccs2
   {
        b1 <- opts$b1      # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2      # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho    # correlation of primary predictor and mediator
        f <- opts$f        # fraction of follow-up times that are uncensored
        fc <- opts$fc      # fraction of failure times that are censored before end of study
        sdx1 <- opts$sdx1  # SD of primary predictor
        sdx2 <- opts$sdx2  # SD of mediator
        ns <- opts$ns      # number of observations in simulated dataset
        S <- matrix(c(sdx1^2, rho*sdx1*sdx2, rho*sdx1*sdx2, sdx2^2), 2, 2)
        X <- matrix(rnorm(2*ns), ns, 2) %*% chol(S)
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.ccs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 25) # (14) cox.bcs2
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1     # prevalence of primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        fc <- opts$fc     # fraction of follow-up times that are censored early
        sdx2 <- opts$sdx2 # SD of mediator
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx2*sqrt(f1/(1-f1))
        mu1 <- rho*sdx2*sqrt((1-f1)/f1)
        sdx2.1 <- sdx2*sqrt(1-rho^2)
        n1 <- round(f1*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rep(0, n0), rnorm(n0, mean=mu0, sd=sdx2.1)),
                   cbind(rep(1, n1), rnorm(n1, mean=mu1, sd=sdx2.1)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bcs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 26) # (14) cox.cbs2
   {
        b1 <- opts$b1     # regression coefficient (log hazard ratio) for primary predictor
        b2 <- opts$b2     # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2     # prevalence of mediator
        rho <- opts$rho   # correlation of primary predictor and mediator
        f <- opts$f       # fraction of follow-up times that are uncensored
        fc <- opts$fc     # fraction of follow-up times that are censored early
        sdx1 <- opts$sdx1 # SD of primary predictor
        ns <- opts$ns     # number of observations in simulated dataset
        mu0 <- -rho*sdx1*sqrt(f2/(1-f2))
        mu1 <- rho*sdx1*sqrt((1-f2)/f2)
        sdx1.2 <- sdx1*sqrt(1-rho^2)
        n1 <- round(f2*ns)
        n0 <- ns-n1
        X <- rbind(cbind(rnorm(n0, mean=mu0, sd=sdx1.2), rep(0, n0)),
                   cbind(rnorm(n1, mean=mu1, sd=sdx1.2), rep(1, n1)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc <- rexp(ns, mean(re)*fc/f)
        } else {
           tc <- rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.cbs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int == 27) # (14) cox.bbs2
   {
        b1 <- opts$b1   # regression coefficient (log hazard ratio) for primary predictor
        f1 <- opts$f1   # prevalence of primary predictor
        b2 <- opts$b2   # regression coefficient (log hazard ratio) for mediator
        f2 <- opts$f2   # prevalence of mediator
        rho <- opts$rho # correlation of primary predictor and mediator
        f <- opts$f     # fraction of follow-up times that are uncensored
        fc <- opts$fc   # fraction of follow-up times that are censored early
        ns <- opts$ns   # number of observations in simulated dataset
        f11 <- rho*sqrt(f1*(1-f1)*f2*(1-f2))+f1*f2
        n11 <- round(ns*f11)
        n10 <- round(ns*(f1-f11))
        n01 <- round(ns*(f2-f11))
        n00 <- ns-n10-n01-n11
        X <- rbind(cbind(rep(0, n00), rep(0, n00)),
                   cbind(rep(1, n10), rep(0, n10)),
                   cbind(rep(0, n01), rep(1, n01)),
                   cbind(rep(1, n11), rep(1, n11)))
        re <- exp(X%*%as.matrix(c(b1, b2)))
        te <- rexp(ns, re)
        if (fc>0) {
           tc = rexp(ns, mean(re)*fc/f)
        } else {
           tc = rep(max(te)+1, ns)
        }
        time <- apply(cbind(te, tc), 1, min)
        event <- ifelse(rank(time)<=round(ns*(f+fc)) & te<tc, 1, 0)
        v2a1 <- survival::coxph(Surv(time, event)~X)$var[2, 2]*ns*f
        desc <- "cox.bbs2"
        d <- (qnorm(alpha)+qnorm(gamma))^2*v2a1/b2^2
        n <- d/f
   }
   if (model.int < 19) list(desc=desc, n=round(n))
   else list(desc=desc, d=round(d), n = round(n))
}

# 23-6-2010 MRC-Epid JHZ
