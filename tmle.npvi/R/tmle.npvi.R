tmle.npvi. <- structure(
    function#Targeted Minimum Loss Estimation of NPVI
### Carries   out  the   targeted  minimum   loss  estimation   (TMLE)   of  a
### non-parametric variable importance measure of a continuous exposure.
    (obs,
### A \code{n  x p} \code{matrix}  or \code{data frame} of  observations, with
### \eqn{p  \ge  3}.  \itemize{  \item{Column  \code{"X"}  corresponds to  the
### continuous  exposure variable  (e.g. DNA  copy  number), or  "cause" in  a
### causal model, with  a reference value \eqn{x_0} equal  to 0.} \item{Column
### \code{"Y"}  corresponds  to the  outcome  variable  (e.g. gene  expression
### level),  or "effect"  in  a  causal model.}  \item{All  other columns  are
### interpreted  as baseline  covariates  \eqn{W} taken  into  account in  the
### definition of the "effect" of \eqn{X} on \eqn{Y}.}}
     f=identity,
### A \code{function} involved in the  definition of the parameter of interest
### \eqn{\psi},  which must  satisfy  \eqn{f(0)=0} (see  Details). Defaults  to
### \code{identity}.
     nMax=30L,
### An \code{integer} (defaults to \code{30L}; \code{10L} is the
### smallest authorized value and we recommend a value less than
### \code{50L} for reasonable computational time) indicating the maximum number of
### observed values of \eqn{X\neq 0} which are used to create the
### supports of the conditional distributions of \eqn{X} given \eqn{W}
### and \eqn{X\neq0} involved in the simulation under \eqn{P_n^k} when
### \code{family} is set to "parsimonious".
     flavor=c("learning", "superLearning"),
### Indicates whether the construction of the relevant features of \eqn{P_n^0}
### and \eqn{P_n^k}, the (non-targeted  yet) initial and (targeted) successive
### updated estimators of the true distribution of \eqn{(W,X,Y)} relies on the
### Super  Learning  methodology   (option  "superLearning")  or  not  (option
### "learning", default  value).  In the former  case, the \code{SuperLearner}
### package is loaded.
     lib=list(),
### A \code{list}  providing the function \code{tmle.npvi}  with the necessary
### algorithms  involved in  the estimation  of the  features of  interest. If
### empty (default) then the default algorithms are used. See \code{Details}.
     nodes=1L,
### An \code{integer},  which indicates how  many nodes should be  involved in
### the  computation  of  the  TMLE  when it  relies  on  the  "superLearning"
### \code{flavor}.  Defaults  to \code{1}.  If larger than  \code{1}, then the
### \code{parallel} package is loaded and a cluster with \code{nodes} nodes is
### created and exploited.
     cvControl=NULL,
### 'NULL' (default value) or an  \code{integer} indicating how many folds are
### involved  in   the  Super   Learning  procedure.   If   \code{flavor}  and
### \code{cvControl} are simultaneously set to "superLearning" and \code{NULL}
### then the Super Learning procedure relies on 10-fold cross-validation.
     family=c("parsimonious", "gaussian"),  #for  'simulateParsimoniouslyXgivenW'   
### Indicates  whether  the  simulation  of the  conditional  distribution  of
### \eqn{X}  given  \eqn{W}  under   \eqn{P_n^k}  (the  initial  estimator  if
### \eqn{k=0} or its  \eqn{k}th update if \eqn{k \ge 1}) should  be based on a
### weighted version  of the  empirical measure (case  "parsimonious", default
### value and faster execution) or on a Gaussian model (case "gaussian").
     cleverCovTheta=FALSE,     #  if TRUE, could be "L2" or "logistic"
### A  \code{logical},  indicating  whether  the  one-step  (if  \code{FALSE},
### default value) or the  two-step (if \code{TRUE}) updating procedure should
### be carried out.
     bound=1,
### A  positive  \code{numeric} (defaults  to  \code{1}),  upper-bound on  the
### absolute value of the fluctuation parameter.
     B=1e5,
### An \code{integer}  (defaults to \code{1e5}) indicating the  sample size of
### the data set  simulated under each \eqn{P_n^k} to  compute an approximated
### value of \eqn{\Psi(P_n^k)}), the parameter of interest at \eqn{P_n^k}. The
### larger \code{B}, the more accurate the approximation.
     trueGMu=NULL,
### Either \code{NULL} (default value) if the \bold{true} conditional
### probability \eqn{g(W)=P(X=0|W)} and conditional expectation
### \eqn{muAux(W)=E_P(X|X\neq 0, W)} are not both known beforehand.
### Otherwise, a \code{list} with tags \code{g} and \code{muAux} and
### respective values the two respective \code{function}s (see the
### related outputs of \code{getSample}).
     iter=5L,
### An \code{integer}  (defaults to \code{5})  indicating a maximal  number of
### updates.
     stoppingCriteria=list(mic=0.01, div=0.01, psi=0.1),
### A   \code{list}    providing   tuning   parameters    for   the   stopping
### criteria.  Defaults  to   \code{list(mic=0.01,  div=0.01,  psi=0.1)},  see
### \code{Details}.
     gmin=5e-2,
### A  positive \code{numeric},  lower-bound on  the  range of  values of  the
### estimated probabilities  \eqn{P_n^k(X=0|W)} that  \eqn{X} be equal  to its
### reference value \code{0} given  \eqn{W}. Defaults to \code{5e-2}, and must
### be smaller than \code{gmax}.
     gmax=.95,
### A positive \code{numeric} smaller  than \code{1}, upper-bound on the range
### of values  of the estimated probabilities  \eqn{P_n^k(X=0|W)} that \eqn{X}
### be  equal to  its  reference  value \code{0}  given  \eqn{W}. Defaults  to
### \code{95e-2}, and must be larger than \code{gmin}.
     mumin=quantile(f(obs[obs[, "X"]!=0, "X"]), type=1, probs=0.01),
### A  \code{numeric}, lower-bound  on the  range of  values of  the estimated
### conditional  expectation  \eqn{E_{P_n^k}(X|X\neq  0,W)} of  \eqn{X}  given
### \eqn{X\neq 0} and \eqn{W}.  Defaults  to the first percentile of \code{X},
### and must be smaller than \code{mumax}.
     mumax=quantile(f(obs[obs[, "X"]!=0, "X"]), type=1, probs=0.99),
### A  \code{numeric}, upper-bound  on the  range of  values of  the estimated
### conditional  expectation  \eqn{E_{P_n^k}(X|X\neq  0,W)} of  \eqn{X}  given
### \eqn{X\neq 0}  and \eqn{W}.  Defaults  to the 99th percentile  of \eqn{X},
### and must be larger than \code{mumin}.
     verbose=FALSE,
### Prescribes the amount of information  output by the function.  Defaults to
### \code{FALSE}.
     tabulate=TRUE, #  If \code{TRUE},  the support of \eqn{P_n^k(X,W) \subset  {(X^i, W^i):  1 \leq  i,j \leq n}}.   Should be  set to  \code{TRUE}. 
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the joint distribution
### of \eqn{(X,W)}  under \eqn{P_n^k} be  a weighted version of  the empirical
### measure (for a faster execution).
     exact=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE}  (its  default  value).   This requires  that  the  successive
### updates  of  the   estimators  of  the  relevant  features   of  the  true
### distribution of \eqn{(W,X,Y)} be derived from the current estimators based
### on exact formulas rather than on their first-order approximations.
     light=TRUE
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
     ) {
      ##alias<< tmle.npvi
      ##seealso<< getSample, getHistory
      
      ##references<< Chambaz, A.,  Neuvial, P., & van der  Laan, M. J. (2012).
      ##Estimation  of  a  non-parametric  variable importance  measure  of  a
      ##continuous exposure. Electronic  journal of statistics, 6, 1059--1099.

      ##references<<  Chambaz, A., Neuvial,  P. (2015).   tmle.npvi: targeted,
      ##integrative search  of associations between  DNA copy number  and gene
      ##expression,   accounting   for   DNA   methylation.   To   appear   in
      ##Bioinformatics Applications Notes.

      ##details<< The  parameter of interest is  defined as \eqn{\psi=\Psi(P)}
      ##with    \deqn{\Psi(P)    =    \frac{E_P[f(X)    *    (\theta(X,W)    -
      ##\theta(0,W))]}{E_P[f(X)^2]},}{\Psi(P)  =  E_P[f(X)  *  (\theta(X,W)  -
      ##\theta(0,W))]  / E_P[f(X)^2],}  with \eqn{P}  the distribution  of the
      ##random vector  \eqn{(W,X,Y)}, \eqn{\theta(X,W) =  E_P[Y|X,W]}, \eqn{0}
      ##\bold{the reference  value for  \eqn{X}}, and \eqn{f}  a user-supplied
      ##function such  that \eqn{f(0)=0} (e.g.,  \eqn{f=identity}, the default
      ##value).
      ##
      ##The TMLE procedure stops when the maximal number of
      ##iterations, \code{iter}, is reached or when at least one of
      ##the following criteria is met:
      ##\itemize{
      ## \item{The empirical mean \eqn{P_n effIC(P_n^{k+1})} of the
      ## efficient influence curve at \eqn{P_n^{k+1}} scaled by the
      ## estimated standard deviation of the efficient influence curve
      ## at \eqn{P_n^{k+1}} is smaller, in absolute value, than
      ## \code{mic}.}
      ## \item{The total variation (TV) distance between P_n^k and
      ## P_n^{k+1} is smaller than \code{div}.}
      ## \item{The change between the successive values
      ## \eqn{Psi(P_n^k)} and \eqn{Psi(P_n^{k+1})} is smaller than
      ## \code{psi}.}
      ##}
      ##
      ##If \code{lib} is an empty list (\code{list()}, default value) then the
      ##default   algorithms   for  the   chosen   \code{flavor}  are   loaded
      ##(\code{learningLib}  when  \code{flavor}   is  set  to  "learning"  or
      ##\code{superLearningLib} when \code{flavor} is set to "superLearning").
      ##A  valid  \code{lib} argument  must  mimick  the  structure of  either
      ##\code{learningLib}    or    \code{superLearningLib},   depending    on
      ##\code{flavor}.
      ##
      ##The  "superLearning"  \code{flavor}  requires the  \code{SuperLearner}
      ##package  and, by  default, the  \code{e1071}, \code{gam},
      ##\code{glmnet}, \code{polspline} and \code{randomForest} packages.
      ## 
      ##If  \code{family}  is set  to  "parsimonious"  (recommended) then  the
      ##packages \code{sgeostat} and \code{geometry} are required.
      
      ## Arguments
      mode <- mode(lib)
      if (mode != "list") {
        throw("Argument 'lib' should be a 'list'");
      }
      test <- setdiff(names(lib), paste0("learn", c("G", "MuAux", "Theta",
                                                    "DevG", "DevMu", "DevTheta",
                                                    "CondExpX2givenW", "CondExpXYgivenW")))
      if (length(test)) {
        throw("Missing element in argument 'lib':", test);
      }
      
      flavor <- match.arg(flavor)
      nMax <- Arguments$getInteger(nMax, c(10, Inf))
      nodes <- Arguments$getInteger(nodes)
      family <- match.arg(family)

      if (B<5e4) {
        warning("Parameter 'B' may be too small; try a larger 'B'")
      }
      
      
      if (length(lib)==0) {
        ## get our library:
        if (flavor=="superLearning") {
          lib <- tmle.npvi::superLearningLib
        } else {
          lib <- tmle.npvi::learningLib
        }
      }
       
      if (flavor=="superLearning") {
        if (is.null(cvControl)) {
          warning("Setting 'V=10' in 'SuperLearner.'")
          cvControl <- SuperLearner::SuperLearner.CV.control(V=10L)
        } else {
          cvControl <- Arguments$getInteger(cvControl,c(2, Inf))
          cvControl <- SuperLearner::SuperLearner.CV.control(V=cvControl)
        }
        if (nodes==1) {
          SuperLearner. <- function(...) {
            SuperLearner::SuperLearner(cvControl=cvControl, ...)
          }
        } else {
          tf <- tempfile("snitch.Rout")
          cl <- parallel::makeCluster(nodes, type="PSOCK", outfile=tf) # can use different types here
          on.exit(parallel::stopCluster(cl))
          ##
          SL.library <- unique(unlist(tmle.npvi::superLearningLib))
          parallel::clusterSetRNGStream(cl, iseed=2343)
          parallel::clusterExport(cl, SL.library)
          SuperLearner. <- function(...) {
            SuperLearner::snowSuperLearner(cluster=cl, cvControl=cvControl, ...)
          }
        }
      } else {
        SuperLearner. <- NULL
        if (nodes>1) {
          warning("Parallel computing not available with 'learning' option")
        }
      }
      
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Declaration
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (any(is.na(obs))) {
        throw("The matrix 'obs' contains at least one 'NA'. This is not allowed.")
      }
      if (is.data.frame(obs)) {
        nms <- colnames(obs)
        varNames <- c("X", "Y")
        m <- match(varNames, nms)
        if (ncol(obs)<3 | any(is.na(m))) {
          throw("The data frame 'obs' must contain at least three columns, including 'X' and 'Y'.")
        }
        XY <- cbind(X=as.numeric(obs[, "X"]), Y=as.numeric(obs[, "Y"]))
        W <- model.matrix(~.-1, obs[, -m])        
        attr(W, "assign") <- NULL
        attr(W, "contrasts") <- NULL
        obs <- cbind(XY, W)
      }

      p0min <- 0.1
      n0min <- p0min*nrow(obs)
      n0 <- sum(obs[, "X"]==0)
      if (n0<n0min) {
        warning("Only ", n0, " out of ", nrow(obs), " observations have 'X==0'. Should 'X' be thresholded?")
      }
      
      npvi <- NPVI(obs=obs, f=f, nMax=nMax, family=family, tabulate=tabulate, 
                   gmin=gmin, gmax=gmax,
                   mumin=mumin, mumax=mumax,
                   thetamin=min(obs[, "Y"]), thetamax=max(obs[, "Y"]),
                   stoppingCriteria=stoppingCriteria)

      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Initialization
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      init(npvi, flavor=flavor,
           learnG=lib$learnG, learnMuAux=lib$learnMuAux, learnTheta=lib$learnTheta,
           bound=bound, B=B,
           light=light,
           trueGMu=trueGMu, 
           SuperLearner.=SuperLearner.,
           verbose=verbose);
      ## rm(learnG);
      conv <- getConv(npvi);
      
      kk <- 1
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Update
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      cat("  iteration ")
      while ((kk <= iter) && (is.na(conv) || !conv)) {
        cat(kk, " ")
        kk <- kk+1

        ## k^th update
        update(npvi, flavor=flavor, learnDevG=lib$learnDevG,
               learnDevMu=lib$learnDevMu, learnDevTheta=lib$learnDevTheta,
               learnCondExpX2givenW=lib$learnCondExpX2givenW,
               learnCondExpXYgivenW=lib$learnCondExpXYgivenW,
               bound=bound, B=B, cleverCovTheta=cleverCovTheta,
               exact=exact, trueGMu=trueGMu,
               SuperLearner.=SuperLearner.,
               verbose=verbose);

        ## check convergence
        updateConv(npvi, B=B)
        conv <- getConv(npvi)
      }
      cat("\n")

      
      ## history <- getHistory(npvi)
      return(npvi)
###   Returns an object of class "NPVI" summarizing the different steps
###   of the  TMLE procedure. The \code{method}  \code{getHistory} outputs the
###   "history" of the procedure  (see \code{getHistory}).  The object notably
###   includes  the following information:
###   \item{obs}{The \code{matrix} of observations  used to carry out the TMLE
###   procedure.  Use the \code{method} \code{getObs} to retrieve it.}
###   \item{psi}{The TMLE of the parameter of interest.  Use the \code{method}
###   \code{getPsi} to retrieve it.}
###   \item{psi.sd}{The  estimated  standard  deviation  of the  TMLE  of  the
###   parameter of interest. Use the \code{method} \code{getPsiSd} to retrieve
###   it.}
    }, ex=function() {
      set.seed(12345)
      ##
      ## Simulating a data set and computing the true value of the parameter
      ##

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
      obs <- sim$obs
      
      ## Adding (dummy) baseline covariates
      V <- matrix(runif(3*nrow(obs)), ncol=3)
      colnames(V) <- paste("V", 1:3, sep="")
      obs <- cbind(V, obs)

      ## Caution! MAKING '0' THE REFERENCE VALUE FOR 'X'
      X0 <- O[2,2]
      obsC <- obs
      obsC[, "X"] <- obsC[, "X"] - X0
      obs <- obsC

      ## True psi and confidence intervals (case 'f=identity')      
      sim <- getSample(1e4, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
      truePsi <- sim$psi

      confInt0 <- truePsi + c(-1, 1)*qnorm(.975)*sqrt(sim$varIC/nrow(sim$obs))
      confInt <- truePsi + c(-1, 1)*qnorm(.975)*sqrt(sim$varIC/nrow(obs))
      cat("\nCase f=identity:\n")
      msg <- paste("\ttrue psi is: ", signif(truePsi, 3), "\n", sep="")
      msg <- paste(msg, "\t95%-confidence interval for the approximation is: ",
                   signif(confInt0, 3), "\n", sep="")
      msg <- paste(msg, "\toptimal 95%-confidence interval is: ",
                   signif(confInt, 3), "\n", sep="")
      cat(msg)

      ##
      ## TMLE procedure
      ##

      ## Running the TMLE procedure
      npvi <- tmle.npvi(obs, f=identity, flavor="learning", B=5e4, nMax=10)

      ## Summarizing its results
      npvi
      setConfLevel(npvi, 0.9)
      npvi
      
      history <- getHistory(npvi)
      print(round(history, 4))

      hp <- history[, "psi"]
      hs <- history[, "sic"]
      hs[1] <- NA
      ics <-  c(-1,1) %*% t(qnorm(0.975)*hs/sqrt(nrow(getObs(npvi))))
      
      pch <- 20
      ylim <- range(c(confInt, hp, ics+hp), na.rm=TRUE)
      
      xs <- (1:length(hs))-1
      plot(xs, hp, ylim=ylim, pch=pch, xlab="Iteration", ylab=expression(psi[n]),
           xaxp=c(0, length(hs)-1, length(hs)-1))
      dummy <- sapply(seq(along=xs), function(x) lines(c(xs[x],xs[x]), hp[x]+ics[, x]))
      
      abline(h=confInt, col=4)
      abline(h=confInt0, col=2)
      
    })

tmle.npvi <- function(obs,         f=identity,        nMax=30L,
                      flavor=c("learning", "superLearning"),  lib=list(),  nodes=1L,  cvControl=NULL,
                      family=c("parsimonious",   "gaussian"),  
                      cleverCovTheta=FALSE, bound=1, B=1e5, trueGMu=NULL,  iter=5L,
                      stoppingCriteria=list(mic=0.01,  div=0.01,  psi=0.1),
                      gmin=5e-2,   gmax=.95,
                      mumin=quantile(f(obs[obs[,  "X"]!=0,   "X"]),  type=1, probs=0.01),
                      mumax=quantile(f(obs[obs[, "X"]!=0, "X"]),  type=1, probs=0.99),
                      verbose=FALSE, tabulate=TRUE, exact=TRUE, light=TRUE) {
  flavor <- match.arg(flavor)
  tmle <- try(tmle.npvi.(obs=obs, f=f, nMax=nMax, flavor=flavor, lib=lib, nodes=nodes, cvControl=cvControl,
                         family=family, cleverCovTheta=cleverCovTheta, bound=bound, B=B,
                         trueGMu=trueGMu, iter=iter,
                         stoppingCriteria=stoppingCriteria, gmin=gmin, gmax=gmax,
                         mumin=mumin, mumax=mumax, verbose=verbose, tabulate=tabulate, exact=exact, light=light))
  failed <- inherits(tmle, "try-error")
  if (flavor=="superLearning" & failed) {
    tmle <- tmle.npvi.(obs=obs, f=f, nMax=nMax, flavor="learning", lib=tmle.npvi::learningLib, nodes=1L,
                       cvControl=cvControl, family=family, cleverCovTheta=cleverCovTheta, bound=bound, B=B, 
                       trueGMu=trueGMu, iter=iter, stoppingCriteria=stoppingCriteria, gmin=gmin, gmax=gmax,
                         mumin=mumin, mumax=mumax, verbose=verbose, tabulate=tabulate, exact=exact, light=light)
    attr(tmle, "flag") <- "Flavor 'superLearning' failed, carried out flavor 'learning' instead."
  }
  return(tmle)
}
## formals(tmle.npvi) <- formals(tmle.npvi.)


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

