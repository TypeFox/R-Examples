getSample <- structure(
    function#Generates Simulated Data
### Generates a  run of  simulated observations of  the form  \eqn{(W,X,Y)} to
### investigate  the  "effect"  of  \eqn{X}  on \eqn{Y}  taking  \eqn{W}  into
### account.
    (n,
### An \code{integer}, the number of observations to be generated.
     O,
### A 3x3 numeric \code{matrix} or \code{data.frame}. Rows are 3 baseline
###   observations used as "class centers" for the simulation.
###   Columns are:
###   \itemize{
###     \item{\eqn{W}, baseline   covariate   (e.g.   DNA   methylation),   or
###     "confounder" in a causal model.}
###     \item{\eqn{X}, continuous exposure variable (e.g. DNA copy number), or
###     "cause" in a causal model ,  with a reference value \eqn{x_0} equal to
###     \code{O[2, "X"]}.}
###     \item{\eqn{Y}, outcome  variable  (e.g.  gene  expression  level),  or
###     "effect" in a causal model.}  }
     lambda0,
### A  \code{function}  that  encodes  the relationship  between  \eqn{W}  and
### \eqn{Y}  in observations  with levels  of \eqn{X}  equal to  the reference
### value \eqn{x_0}.
     p=rep(1/3, 3),
### A \code{vector} of  length 3 whose entries sum to  1. Entry \code{p[k]} is
### the probability that each observation  belong to class \code{k} with class
### center \code{O[k, ]}. Defaults to \code{rep(1/3, 3)}.
     omega=rep(1/2, 3),
### A \code{vector}  of length 3 with positive  entries. Entry \code{omega[k]}
### is  the  standard deviation  of  \eqn{W} in  class  \code{k}  (on a  logit
### scale). Defaults to \code{rep(1/2, 3)}.
     Sigma1=matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), 2, 2),
### A 2x2 covariance  \code{matrix} of the random vector  \eqn{(X,Y)} in class
### \code{k=1}, assumed  to be bivariate Gaussian with  mean \code{O[1, c("X",
### "Y")]}. Defaults to \code{matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), 2, 2)}.
     sigma2=omega[2]/5,
### A positive \code{numeric}, the variance  of the random variable \eqn{Y} in
### class \code{k=2},  assumed to be univariate Gaussian  with mean \code{O[2,
### "Y"]}. Defaults to \code{omega[2]/5}.
     Sigma3=Sigma1,
### A 2x2 covariance  \code{matrix} of the random vector  \eqn{(X,Y)} in class
### \code{k=3}, assumed  to be bivariate Gaussian with  mean \code{O[3, c("X",
### "Y")]}. Defaults to \code{Sigma1}.
     f=identity,
### A \code{function} involved in the  definition of the parameter of interest
### \eqn{\psi},  which must  satisfy  \eqn{f(0)=0} (see  Details). Defaults  to
### \code{identity}.
     verbose=FALSE
### Prescribes the amount of information  output by the function.  Defaults to
### \code{FALSE}.
     ) {
      ##references<< Chambaz, A., Neuvial, P., & van der Laan, M. J. (2012).
      ##Estimation of a non-parametric variable importance measure of a
      ##continuous exposure. Electronic journal of statistics, 6, 1059--1099.
      ##details<< The  parameter of interest is  defined as \eqn{\psi=\Psi(P)}
      ##   with   \deqn{\Psi(P)   =   \frac{E_P[f(X-x_0)  *   (\theta(X,W)   -
      ##   \theta(x_0,W))]}{E_P[f(X-x_0)^2]},}{\Psi(P)    =   E_P[f(X-x_0)   *
      ##   (\theta(X,W) - \theta(x_0,W))] / E_P[f(X-x_0)^2],} with \eqn{P} the
      ##   distribution of the random vector \eqn{(W,X,Y)}, \eqn{\theta(X,W) =
      ##   E_P[Y|X,W]}, \eqn{x_0} the reference value for \eqn{X}, and \eqn{f}
      ##   a   user-supplied   function    such   that   \eqn{f(0)=0}   (e.g.,
      ##   \eqn{f=identity},  the  default value).   The  value \eqn{\psi}  is
      ##   obtained  using   the  \bold{known}  \eqn{\theta}   and  the  joint
      ##   empirical  distribution of  \eqn{(X,W)} based  on the  same  run of
      ##   observations  as  in  \code{obs}.   Seeing  \eqn{W, X,  Y}  as  DNA
      ##   methylation, DNA copy number and gene expression, respectively, the
      ##   simulation scheme implements the following constraints:
      ##   \itemize{
      ##     \item There are two or  three copy number classes: normal regions
      ##       (\code{k=2}), and  regions of  copy number gains  and/or losses
      ##       (\code{k=1} and/or \code{k=3}).
      ##     \item In  normal regions,  gene expression levels  are negatively
      ##       correlated with DNA methylation.
      ##     \item  In regions  of  copy number  alteration,  copy number  and
      ##     expression are positively correlated.  }
      

      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Validate arguments
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Argument 'n':
      n <- Arguments$getNumeric(n);
      
      ## Argument 'O':
      if (!is.matrix(O) && !is.data.frame(O)) {
        throw("Argument 'O' should be a matrix or a data.frame: ", mode(O)[1]);
      }
      varNames <- c("W", "X", "Y");
      nms <- colnames(O);
      m <- match(varNames, nms);
      idxs <- which(is.na(m));
      if (length(idxs)) {
        throw("Missing observation:", varNames[idxs]);
      }
      O <- O[, varNames];
      
      ## Argument 'lambda0':
      mode <- mode(lambda0);
      if (mode != "function") {
        throw("Argument 'lambda0' should be of mode 'function', not '", mode);
      }

      ## Argument 'p':
      p <- Arguments$getNumerics(p, range=c(0, 1));
      if (sum(p) != 1) {
        throw("Elements of 'p' should sum to 1");
      }

      ## Argument 'omega':
      omega <- Arguments$getNumerics(omega, range=c(0, Inf));

      ## Argument 'Sigma1':
      Sigma1 <- Arguments$getNumerics(Sigma1);
      if ((nrow(Sigma1) != 2) || (ncol(Sigma1) != 2)) {
        throw("Argument 'Sigma1' should be a 2x2 matrix");
      }

      ## Argument 'sigma2':
      sigma2 <- Arguments$getNumeric(sigma2);
      if (is.na(sigma2)) {
        sigma2 <- mean(X^2);
      }

      ## Argument 'Sigma3':
      Sigma3 <- Arguments$getNumerics(Sigma3);
      if ((nrow(Sigma3) != 2) || (ncol(Sigma3) != 2)) {
        throw("Argument 'Sigma3' should be a 2x2 matrix");
      }

      ## Argument 'f':
      if (!((mode(f)=="function") && (f(0)==0))) {
        throw("Argument 'f' must be a function such that f(0)=0.")
      }

      ## Argument 'verbose':
      verbose <- Arguments$getVerbose(verbose);
      verbose <- less(verbose, 10);

      U <- findInterval(runif(n), cumsum(c(0, p)));
      verbose && print(verbose, table(U))
      ##  W <- O[U, "W"] + rnorm(n, mean=rep(0, n), sd=omega[U]); ## doesn't work: W should be in [0,1]
      logit <- qlogis
      expit <- plogis
      
      W <- expit(logit(O[U, "W"]) + rnorm(n, mean=rep(0, n), sd=omega[U]));
      X <- rep(NA, n);
      Y <- rep(NA, n);

      verbose && enter(verbose, "Simulated copy number and expression data for class 1");
      idx <- which(U == 1);
      verbose && str(verbose, idx);
      if (length(idx)) {
        mu <- as.matrix(O[1, c("X", "Y"), drop=FALSE]);
        verbose && cat(verbose, "mu:");
        verbose && print(verbose, mu);
        XY <- mvrnorm(length(idx), mu=mu, Sigma=Sigma1);
        X[idx] <- XY[, 1];
        Y[idx] <- XY[, 2];
      }
      verbose && exit(verbose);
      
      verbose && enter(verbose, "Simulated copy number and expression data for class 2");
      idx <- which(U == 2);
      verbose && str(verbose, idx);
      if (length(idx)) {
        X[idx] <- O[2, "X"];
        Y[idx] <- O[2, "Y"] + lambda0(W[idx]) +
            rnorm(length(idx), mean=0, sd=sigma2);
      }
      verbose && exit(verbose);
      
      verbose && enter(verbose, "Simulated copy number and expression data for class 3");
      idx <- which(U == 3);
      verbose && str(verbose, idx);
      if (length(idx)) {
        mu <- as.matrix(O[3, c("X", "Y"), drop=FALSE]);
        verbose && cat(verbose, "mu:");
        verbose && print(verbose, mu);
        XY <- mvrnorm(length(idx), mu=mu, Sigma=Sigma3);
        X[idx] <- XY[, 1];
        Y[idx] <- XY[, 2];
      }
      verbose && exit(verbose);
      

      condProbUW <- function(W) {
        u <- 1
        cpu1 <- rep(0, length(W))
        if (p[u]!=0) {
          cpu1 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
        }

        u <- 2
        cpu2 <- rep(0, length(W))
        if (p[u]!=0) {
          cpu2 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
        }

        u <- 3
        cpu3 <- rep(0, length(W))
        if (p[u]!=0) {
          cpu3 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) / omega[u]
        }
        out <- cbind(cpu1, cpu2, cpu3)
        ## make sure that they sum to 1
        Z <- apply(out, 1, sum, na.rm = TRUE)
        out <- out/Z
        out
      }

      g <- function(W) {
        condProbUW(W)[, 2]
      }
      
      condProbUXW <- function(XW) {
        X <- XW[, 1]
        W <- XW[, 2]
        
        u <- 1
        cpu1 <- rep(0, nrow(XW))
        if (p[u]!=0) {
          cpu1 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
                  dnorm(X, mean=O[u, "X"], sd=Sigma1[1, 1]) / omega[u]
        }
        
        u <- 2
        cpu2 <- rep(0, nrow(XW))
        if (p[u]!=0) {
          cpu2 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
                  (X==O[u, "X"]) / omega[u]
        }

        u <- 3
        cpu3 <- rep(0, nrow(XW))
        if (p[u]!=0) {
          cpu3 <- p[u] *
              dnorm(logit(W), mean=logit(O[u, "W"]), sd=omega[u])/(W*(1-W)) *
                  dnorm(X, mean=O[u, "X"], sd=Sigma3[1, 1]) / omega[u]
        }
        out <- cbind(cpu1, cpu2, cpu3)
        ## make sure that they sum to 1
        Z <- apply(out, 1, sum, na.rm = TRUE)
        out <- out/Z
        out
      }
      
      ## Conditional expectation of X given W
      mu <- function(W) {
        vec <- f(O[, "X"] - O[2, "X"])
        vec[is.na(vec)] <- 0 ## any real number would do here!
        res <- condProbUW(W) %*% vec
        dim(res) <- NULL
        res
      }

      muAux <- function(W) {
        mu(W)/(1-g(W))
      }
      
      ## Conditional expectation of Y given (W,X)
      theta <- function(XW) {
        XW[, "X"] <- XW[, "X"] + O[2,"X"]
        X <- XW[, 1]
        W <- XW[, 2]
        
        dummy <- 1
        param <- c(Sigma1[1,2]/Sigma1[1,1],
                   dummy,
                   Sigma3[1,2]/Sigma3[1,1])
        ## if (X_1,X_2)~N(mu, S) then E(X_2|X_1)=mu_2+S_12/S_11*(X_1-mu_1)
        cpu <- condProbUXW(XW)
        
        res <- 0
        res1 <- (O[1, "Y"] + param[1] * (X-O[1, "X"])) * cpu[, 1]
        if (p[1] != 0) {
          res <- res + res1
        }
        
        res2 <- (O[2, "Y"] + lambda0(W)) * cpu[, 2]
        if (p[2] != 0) {
          res <- res+res2
        }
        
        res3 <- (O[3, "Y"] + param[3] * (X-O[3, "X"])) * cpu[, 3]
        if (p[3] != 0) {
          res <- res + res3
        }
        
        
        dim(res) <- NULL
        res
      }
      
      ## function 'theta0'
      theta0 <- function(W) {
        XW <- cbind(X=0, W=W);
        theta(XW);
      }

      if (identical(f, identity)) {
        sigma2 <- p * c(Sigma1[1, 1] + (O[1, "X"]-O[2, "X"])^2,
                        0,
                        Sigma3[1, 1] + (O[3, "X"]-O[2, "X"])^2)
        sigma2 <- sum(sigma2, na.rm=TRUE)
      } else {
        foo <- function(x, ...) {
          f(x)^2*dnorm(x, ...)
        }
        m <- O[1, "X"]-O[2, "X"]
        sigma21 <- NA
        if (!is.na(m)) {
          sigma21 <- integrate(foo, lower=-Inf, upper=+Inf,
                               mean=m, sd=sqrt(Sigma1[1, 1]))$value
        } 
        m <- O[3, "X"]-O[2, "X"]
        sigma23 <- NA
        if (!is.na(m)) {
          sigma23 <- integrate(foo, lower=-Inf, upper=+Inf,
                               mean=m, sd=sqrt(Sigma3[1, 1]))$value
        }
        sigma2 <- sum(p*c(sigma21, 0, sigma23), na.rm=TRUE)
      }
      
      theX0 <- O[2,"X"];
      obs <- cbind(W,X,Y)
      rownames(obs) <- NULL
      rm(W,X,Y,U)

      obsC <- obs;
      obsC[, "X"] <- obsC[, "X"] - theX0;
      
      fX <- function(obs) {
        f(obs[, "X"]);
      }
      
      res <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obsC, sigma2=sigma2);
      truePsi <- res$mean
      sd.truePsi <- res$sd
      
      effIC <- function(obs) {
        ## Argument 'obs':
        if (FALSE) {
          if (!is.matrix(obs)) {
            throw("Argument 'obs' should be a matrix: ", mode(obs)[1]);
          }
        }
        obs <- validateArgumentObs(obs, allowIntegers=FALSE);
        
        W <- obs[, "W"]
        verbose && str(verbose, W);
        X <- f(obs[, "X"])
        ## CAUTION! 'X' is 'f(X)' indeed and 'theta' does apply to 'obs[, c("X", "W")]'...
        Y <- obs[, "Y"]
        eic1 <- X*(theta(obs[, c("X", "W")]) - theta(cbind(X=0, W=W)) - X*truePsi)
        eic2 <- (Y - theta(obs[, c("X", "W")])) * (X - mu(W)/g(W)*(X==0))
        out <- 1/sigma2*(eic1+eic2)
        dim(out) <- NULL
        out
      }

      varIC <- var(effIC(obsC));
      
      rm(obsC);
      
      res <- list(
          obs=obs,
          psi=truePsi,
          g=g,
          mu=mu,
          muAux=muAux,
          theta=theta,
          theta0=theta0,
          sigma2=sigma2,
          effIC=effIC,
          varIC=varIC);
      return(res);
###   Returns a list with the following tags:
###   \item{obs}{A \code{matrix} of \code{n} observations. The \code{c(W,X,Y)}
###   colums of \code{obs} have the  same interpretation as the columns of the
###   input argument \code{O}.}
###   \item{psi}{A  \code{numeric}, approximated value  of the  true parameter
###   \eqn{\psi} obtained  using the  \bold{known} \eqn{\theta} and  the joint
###   empirical  distribution  of  \eqn{(X,W)}   based  on  the  same  run  of
###   observations  as in  \code{obs}. The  larger the  sample size,  the more
###   accurate the approximation.}
###   \item{g}{A   \code{function},  the  \bold{known}   positive  conditional
###   probability   \eqn{P(X=x_0|W)}.}
###   \item{mu}{A  \code{function}, the  \bold{known}  conditional expectation
###   \eqn{E_P(X|W)}.}
###   \item{muAux}{A \code{function}, the \bold{known} conditional expectation
###   \eqn{E_P(X|X\neq x_0, W)}.}
###   \item{theta}{A \code{function}, the \bold{known} conditional expectation
###   \eqn{E_P(Y|X,W)}.}
###   \item{theta0}{A    \code{function},    the   \bold{known}    conditional
###   expectation \eqn{E_P(Y|X=x_0,W)}.}
###   \item{sigma2}{A  positive \code{numeric},  the  \bold{known} expectation
###   \eqn{E_P(f(X-x_0)^2)}.}
###   \item{effIC}{A  \code{function},  the  \bold{known} efficient  influence
###   curve of  the functional \eqn{\Psi} at \eqn{P},  \bold{assuming} that the
###   reference value \eqn{x_0=0}.}
###   \item{varIC}{A  positive  \code{numeric},   approximated  value  of  the
###   variance of the efficient influence curve of the functional \eqn{\Psi} at
###   \eqn{P}  and  evaluated at  \eqn{O},  obtained  using  the same  run  of
###   observations  as in  \code{obs}. The  larger the  sample size,  the more
###   accurate the approximation. }
    }, ex=function() {
      ## Parameters for the simulation (case 'f=identity')
      O <- cbind(W=c(0.05218652, 0.01113460),
                 X=c(2.722713, 9.362432),
                 Y=c(-0.4569579, 1.2470822))
      O <- rbind(NA, O)
      lambda0 <- function(W) {-W}
      p <- c(0, 1/2, 1/2)
      omega <- c(0, 3, 3)
      S <- matrix(c(10, 1, 1, 0.5), 2 ,2)

      ## Simulating a data set of 200 i.i.d. observations
      sim <- getSample(2e2, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
      str(sim)
      
      obs <- sim$obs
      head(obs)
      pairs(obs)

      ## Adding (dummy) baseline covariates
      V <- matrix(runif(3*nrow(obs)), ncol=3)
      colnames(V) <- paste("V", 1:3, sep="")
      obsV <- cbind(V, obs)
      head(obsV)

      ## True psi and confidence intervals (case 'f=identity')      
      sim01 <- getSample(1e4, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
      truePsi1 <- sim01$psi

      confInt01 <- truePsi1+c(-1, 1)*qnorm(.975)*sqrt(sim01$varIC/nrow(sim01$obs))
      confInt1 <- truePsi1+c(-1, 1)*qnorm(.975)*sqrt(sim01$varIC/nrow(obs))
      msg <- "\nCase f=identity:\n"
      msg <- c(msg, "\ttrue psi is: ", paste(signif(truePsi1, 3)), "\n")
      msg <- c(msg, "\t95%-confidence interval for the approximation is: ", 
               paste(signif(confInt01, 3)), "\n")
      msg <- c(msg, "\toptimal 95%-confidence interval is: ", 
               paste(signif(confInt1, 3)), "\n")
      cat(msg)

      ## True psi and confidence intervals (case 'f=atan')
      f2 <- function(x) {1*atan(x/1)}
      sim02 <- getSample(1e4, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S, f=f2);
      truePsi2 <- sim02$psi;
      
      confInt02 <- truePsi2+c(-1, 1)*qnorm(.975)*sqrt(sim02$varIC/nrow(sim02$obs))
      confInt2 <- truePsi2+c(-1, 1)*qnorm(.975)*sqrt(sim02$varIC/nrow(obs))

      msg <- "\nCase f=atan:\n"
      msg <- c(msg, "\ttrue psi is: ", paste(signif(truePsi2, 3)), "\n")
      msg <- c(msg, "\t95%-confidence interval for the approximation is: ", 
               paste(signif(confInt02, 3)), "\n")
      msg <- c(msg, "\toptimal 95%-confidence interval is: ", 
               paste(signif(confInt2, 3)), "\n")
      cat(msg)
    })

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

