

#'@rdname mc.est
#'@method mc.est CBData
#'@export

mc.est.CBData <- function(object, ...){
  cbdata <- object[object$Freq>0, ]
  #by trt
  do.est.fun <- function(x){
    est <- .Call("ReprodEstimates", as.integer(x$ClusterSize), as.integer(x$NResp), 
                             as.integer(x$Freq),PACKAGE="CorrBin")
    est <- cbind(c(1,rep(NA,nrow(est)-1)), est) 
    idx <- upper.tri(est,diag=TRUE)
    est.d <- data.frame(Prob=est[idx], ClusterSize=as.integer(col(est)[idx]-1), 
                        NResp=as.integer(row(est)[idx]-1),
                        Trt=x$Trt[1])
    est.d}  
  
  est.list <- by(cbdata, list(Trt=cbdata$Trt), do.est.fun)
  do.call(rbind, est.list)}


#'@rdname mc.test.chisq
#'@method mc.test.chisq CBData
#'@export

mc.test.chisq.CBData <- function(object,...){
  cbdata <- object[object$Freq>0, ]
 
  get.T <- function(x){
      max.size <- max(x$ClusterSize)
      scores <- (1:max.size) - (max.size+1)/2
      p.hat <- with(x, sum(Freq*NResp) / sum(Freq*ClusterSize))
      rho.hat <- with(x, 1-sum(Freq*(ClusterSize-NResp)*NResp/ClusterSize) / 
          (sum(Freq*(ClusterSize-1))*p.hat*(1-p.hat)))  #Fleiss-Cuzick estimate
      c.bar <- with(x, sum(Freq*scores[ClusterSize]*ClusterSize) / sum(Freq*ClusterSize))
      T.center <- with(x, sum(Freq*(scores[ClusterSize]-c.bar)*NResp))
      Var.T.stat <-  with(x, 
         p.hat*(1-p.hat)*sum(Freq*(scores[ClusterSize]-c.bar)^2*ClusterSize*(1+(ClusterSize-1)*rho.hat)))
      X.stat <- (T.center)^2/Var.T.stat
      X.stat}
      
   chis <- by(cbdata, cbdata$Trt, get.T)
   chis <- chis[1:length(chis)]
   chi.list <- list(chi.sq=chis, p=pchisq(chis, df=1, lower.tail=FALSE))
   overall.chi <- sum(chis)
   overall.df <- length(chis)
   list(overall.chi=overall.chi, overall.p=pchisq(overall.chi, df=overall.df, lower.tail=FALSE), 
        individual=chi.list)
}

#'Order-restricted MLE assuming marginal compatibility
#'
#'\code{SO.mc.est} computes the nonparametric maximum likelihood estimate of
#'the distribution of the number of responses in a cluster \eqn{P(R=r|n)} under
#'a stochastic ordering constraint. Umbrella ordering can be specified using
#'the \code{turn} parameter.
#'
#'Two different algorithms: EM and ISDM are implemented. In general, ISDM (the
#'default) should be faster, though its performance depends on the tuning
#'parameter \code{max.directions}: values that are too low or too high slow the
#'algorithm down.
#'
#'\code{SO.mc.est} allows extension to an umbrella ordering: \eqn{D_1 \geq^{st}
#'\cdots \geq^{st} D_k \leq^{st} \cdots \leq^{st} D_n}{D_1 >= \ldots >= D_k <=
#'\ldots <= D_n} by specifying the value of \eqn{k} as the \code{turn}
#'parameter. This is an experimental feature, and at this point none of the
#'other functions can handle umbrella orderings.
#'
#'@useDynLib CorrBin
#'@export
#'@param cbdata an object of class \code{\link{CBData}}.
#'@param turn integer specifying the peak of the umbrella ordering (see
#'Details). The default corresponds to a non-decreasing order.
#'@param control an optional list of control settings, usually a call to
#'\code{\link{soControl}}.  See there for the names of the settable control
#'values and their effect.
#'@return A list with components:
#'
#'Components \code{Q} and \code{D} are unlikely to be needed by the user.
#'@return \item{MLest}{data frame with the maximum likelihood estimates of
#'\eqn{P(R_i=r|n)}}
#'@return \item{Q}{numeric matrix; estimated weights for the mixing distribution}
#'@return \item{D}{numeric matrix; directional derivative of the log-likelihood}
#'@return \item{loglik}{the achieved value of the log-likelihood}
#'@return \item{converge}{a 2-element vector with the achived relative error and
#'the performed number of iterations}
#'@author Aniko Szabo, aszabo@@mcw.edu
#'@seealso \code{\link{soControl}}
#'@references Szabo A, George EO. (2009) On the Use of Stochastic Ordering to
#'Test for Trend with Clustered Binary Data. \emph{Biometrika}
#'@keywords nonparametric models
#'@examples
#'
#'  data(shelltox)
#'  ml <- SO.mc.est(shelltox, control=soControl(eps=0.01, method="ISDM"))
#'  attr(ml, "converge")
#'  
#'  require(lattice)
#'  panel.cumsum <- function(x,y,...){
#'    x.ord <- order(x)
#'    panel.xyplot(x[x.ord], cumsum(y[x.ord]), ...)}
#'
#'  xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=ml, type="s",
#'       panel=panel.superpose, panel.groups=panel.cumsum,
#'       as.table=TRUE, auto.key=list(columns=4, lines=TRUE, points=FALSE),
#'       xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
#'       ylim=c(0,1.1), main="Stochastically ordered estimates\n with marginal compatibility")
#'


SO.mc.est <- function(cbdata, turn=1, control=soControl()){ 
  tab <- xtabs(Freq~factor(ClusterSize,levels=1:max(ClusterSize))+
                factor(NResp,levels=0:max(ClusterSize))+Trt, data=cbdata)
        size <- dim(tab)[1]
        ntrt <- dim(tab)[3]
        ntot <- sum(tab)
  storage.mode(tab) <- "double"
        Q <- array(0, dim=rep(size+1,ntrt))
        storage.mode(Q) <- "double"
          
        S <- DownUpMatrix(size, ntrt, turn)
        storage.mode(S) <- "integer"
        
        if ((control$start=="H0")&(control$method=="EM")){
          warning("The EM algorithm can only use 'start=uniform'. Switching options.")
          start <- "uniform"
  }
        if (control$start=="H0"){
          const.row <- matrix(0:size, nrow=size+1, ncol=ntrt)
          Q[const.row+1] <- 1/(size+1)
                }
        else {  #start=="uniform"
          Q[S+1] <- 1/(nrow(S))
    }
         
  res0 <- switch(control$method,
      EM = .Call("MixReprodQ", Q, S, tab, as.integer(control$max.iter), as.double(control$eps), 
                                   as.integer(control$verbose), PACKAGE="CorrBin"),
      ISDM = .Call("ReprodISDM", Q, S, tab, as.integer(control$max.iter), as.integer(control$max.directions),
                   as.double(control$eps),  as.integer(control$verbose), PACKAGE="CorrBin"))
 
  names(res0) <- c("MLest","Q","D","loglik", "converge")
  names(res0$converge) <- c("rel.error", "n.iter")
  res <- res0$MLest
 
  dimnames(res) <- list(NResp=0:size, ClusterSize=1:size, Trt=1:ntrt)
  res <- as.data.frame.table(res)
  names(res) <- c("NResp","ClusterSize","Trt","Prob") 
  res$NResp  <- as.numeric(as.character(res$NResp))
  res$ClusterSize  <- as.numeric(as.character(res$ClusterSize))
  res <- res[res$NResp <= res$ClusterSize,]
  levels(res$Trt) <- levels(cbdata$Trt)
  
  attr(res, "loglik") <- res0$loglik
  attr(res, "converge") <- res0$converge
  res
}

#'Control values for order-constrained fit
#'
#'The values supplied in the function call replace the defaults and a list with
#'all possible arguments is returned.  The returned list is used as the control
#'argument to the \code{\link{mc.est}}, \code{\link{SO.LRT}}, and
#'\code{\link{SO.trend.test}} functions.
#'
#'@export
#'@param method a string specifying the maximization method
#'@param eps a numeric value giving the maximum absolute error in the
#'log-likelihood
#'@param max.iter an interger specifying the maximal number of iterations
#'@param max.directions an integer giving the maximal number of directions
#'considered at one step of the ISDM method.  If zero or negative, it is set to
#'the number of non-empty cells. A value of 1 corresponds to the VDM algorithm.
#'@param start a string specifying the starting setup of the mixing
#'distribution; "H0" puts weight only on constant vectors (corresponding to the
#'null hypothesis of no change), "uniform" puts equal weight on all elements.
#'Only a "uniform" start can be used for the "EM" algorithm.
#'@param verbose a logical value; if TRUE details of the optimization are
#'shown.
#'@return a list with components for each of the possible arguments.
#'@author Aniko Szabo aszabo@@mcw.edu
#'@seealso \code{\link{mc.est}}, \code{\link{SO.LRT}},
#'\code{\link{SO.trend.test}}
#'@keywords models
#'@examples
#'
#'# decrease the maximum number iterations and
#'# request the "EM" algorithm
#' soControl(method="EM", max.iter=100)
#'

soControl <- function(method=c("ISDM","EM"), eps=0.005, max.iter=5000, 
      max.directions=0, start=ifelse(method=="ISDM", "H0", "uniform"), verbose=FALSE){
  method <- match.arg(method)
  start <- match.arg(start, c("uniform","H0"))
  list(method = match.arg(method), eps = eps, max.iter = max.iter,
       max.directions = max.directions, start=start, verbose = verbose)
}

#'Internal CorrBin objects
#'
#'Internal CorrBin objects.
#'
#'These are not to be called by the user.
#'
#'@rdname CorrBin-internal
#'@aliases .required DownUpMatrix
#'@keywords internal

DownUpMatrix <- function(size, ntrt, turn){
  if ((turn<1)|(turn>ntrt)) stop("turn should be between 1 and ntrt")
  
    if (turn==1){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
       return(res)
    }
    if (turn==ntrt){
       res <- .Call("makeSmatrix", as.integer(size), as.integer(ntrt),PACKAGE="CorrBin")
       return(size - res)
    }
  
  
    res1 <- .Call("makeSmatrix", as.integer(size), as.integer(turn),PACKAGE="CorrBin")
    res1 <- size - res1;
  
  res2list <- list()
  for (sq in 0:size){
    
      S <- .Call("makeSmatrix", as.integer(size-sq), as.integer(ntrt-turn),PACKAGE="CorrBin")
      res2list <- c(res2list, list(sq+S))
    
  }
  
    res1list <- by(res1, res1[,turn], function(x)x)
    res <- mapply(merge, res1list, res2list, MoreArgs=list(by=NULL), SIMPLIFY=FALSE)
    res <- data.matrix(do.call(rbind, res))
    rownames(res) <- NULL
    colnames(res) <- NULL
    res
  
}

#'Likelihood-ratio test statistic
#'
#'\code{SO.LRT} computes the likelihood ratio test statistic for stochastic
#'ordering against equality assuming marginal compatibility for both
#'alternatives. Note that this statistic does not have a
#'\eqn{\chi^2}{chi-squared} distribution, so the p-value computation is not
#'straightforward. The \code{\link{SO.trend.test}} function implements a
#'permutation-based evaluation of the p-value for the likelihood-ratio test.
#'
#'@export
#'@param cbdata a \code{CBData} object
#'@param control an optional list of control settings, usually a call to
#'\code{\link{soControl}}.  See there for the names of the settable control
#'values and their effect.
#'@return The value of the likelihood ratio test statistic is returned with two
#'attributes:
#'@return \item{ll0}{the log-likelihood under \eqn{H_0}{H0} (equality)}
#'@return \item{ll1}{the log-likelihood under \eqn{H_a}{Ha} (stochastic order)}
#'@author Aniko Szabo
#'@seealso \code{\link{SO.trend.test}}, \code{\link{soControl}}
#'@keywords htest nonparametric
#'@examples
#'
#'data(shelltox)
#'LRT <- SO.LRT(shelltox, control=soControl(max.iter = 100, max.directions = 50))
#'LRT
#'

SO.LRT <- function(cbdata, control=soControl()){
        # LL under null hypothesis of equality (+ reproducibility)
        a <- with(cbdata, aggregate(Freq, list(ClusterSize=ClusterSize,NResp=NResp), sum))
        names(a)[names(a)=="x"] <- "Freq"
        a$ClusterSize <- as.integer(as.character(a$ClusterSize))
        a$NResp <- as.integer(as.character(a$NResp))
        a$Trt <- 1
    class(a) <- c("CBData", "data.frame")
                       
        b <- mc.est(a)
  b <- merge(cbdata, b, all.x=TRUE, by=c("ClusterSize","NResp"))
  ll0 <- with(b, sum(Freq*log(Prob)))
  
  # LL under alternative hypothesis of stoch ordering (+ reproducibility)
  res <- SO.mc.est(cbdata, control=control)
  ll1 <- attr(res, "loglik")
  lrt <- 2*(ll1 - ll0)
  attr(lrt, "ll0") <- ll0
  attr(lrt, "ll1") <- ll1
  lrt
 }
  

#'Likelihood ratio test of stochastic ordering
#'
#'Performs a likelihood ratio test of stochastic ordering versus equality using
#'permutations to estimate the null-distribution and the p-value.  If only the
#'value of the test statistic is needed, use \code{\link{SO.LRT}} instead.
#'
#'The test is valid only under the assumption that the cluster-size
#'distribution does not depend on group. During the estimation of the
#'null-distribution the group assignments of the clusters are permuted keeping
#'the group sizes constant; the within-group distribution of the cluster-sizes
#'will vary randomly during the permutation test.
#'
#'The default value of \code{R} is probably too low for the final data
#'analysis, and should be increased.
#'
#'@import boot
#'@export
#'@param cbdata a \code{\link{CBData}} object.
#'@param R an integer -- the number of random permutations for estimating the
#'null distribution.
#'@param control an optional list of control settings, usually a call to
#'\code{\link{soControl}}.  See there for the names of the settable control
#'values and their effect.
#'@return A list with the following components
#'@return \item{LRT}{the value of the likelihood ratio test statistic. It has two
#'attributes: \code{ll0} and \code{ll1} - the values of the log-likelihood
#'under \eqn{H_0}{H0} and \eqn{H_a}{Ha} respectively.}
#'@return \item{p.val}{the estimated one-sided p-value.}
#'@return \item{boot.res}{an object of class "boot" with the detailed results of
#'the permutations.  See \code{\link[boot]{boot}} for details.}
#'@author Aniko Szabo, aszabo@@mcw.edu
#'@seealso \code{\link{SO.LRT}} for calculating only the test statistic,
#'\code{\link{soControl}}
#'@references Szabo A, George EO. (2009) On the Use of Stochastic Ordering to
#'Test for Trend with Clustered Binary Data.
#'@keywords htest nonparametric
#'@examples
#'
#'data(shelltox)
#'set.seed(45742)
#'sh.test <- SO.trend.test(shelltox, R=10, control=soControl(eps=0.1, max.directions=25)) 
#'sh.test
#'
#'#a plot of the resampled LRT values
#'#would look better with a reasonable value of R
#' null.vals <- sh.test$boot.res$t[,1]
#' hist(null.vals, breaks=10,  freq=FALSE, xlab="Test statistic", ylab="Density", 
#'      main="Simulated null-distribution", xlim=range(c(0,20,null.vals)))
#' points(sh.test$LRT, 0, pch="*",col="red", cex=3)
#'

SO.trend.test <- function(cbdata, R=100, control=soControl()){
        dat2 <- cbdata[rep(1:nrow(cbdata), cbdata$Freq),]  #each row is one sample
        dat2$Freq <- NULL
        
        boot.LRT.fun <- function(dat, idx){
          dat.new <- cbind(dat[idx, c("ClusterSize","NResp")], Trt=dat$Trt)   #rearrange clusters
                dat.f <- aggregate(dat.new$Trt, 
                          list(Trt=dat.new$Trt, ClusterSize=dat.new$ClusterSize, NResp=dat.new$NResp), length)
          names(dat.f)[names(dat.f)=="x"] <- "Freq"
    dat.f$ClusterSize <- as.numeric(as.character(dat.f$ClusterSize))
          dat.f$NResp <- as.numeric(as.character(dat.f$NResp))
    class(dat.f) <- c("CBData", class(dat.f))
                    
    stat <- SO.LRT(dat.f, control=control)
    stat}        
        
   res <- boot(dat2, boot.LRT.fun, R=R, sim="permutation")
     
   p <- mean(res$t[,1] >= res$t0)
   LRT <- res$t0
   list(LRT=LRT, p.val=p, boot.res=res)}           

#'Test for increasing trend with correlated binary data
#'
#'The \code{trend.test} function provides a common interface to the trend tests
#'implemented in this package: \code{\link{SO.trend.test}},
#'\code{\link{RS.trend.test}}, and \code{\link{GEE.trend.test}}. The details of
#'each test can be found on their help page.
#'
#'@export
#'@param cbdata a \code{\link{CBData}} object
#'@param test character string defining the desired test statistic. "RS"
#'performs the Rao-Scott test (\code{\link{RS.trend.test}}), "SO" performs the
#'stochastic ordering test (\code{\link{SO.trend.test}}), "GEE", "GEEtrend",
#'"GEEall" perform the GEE-based test (\code{\link{GEE.trend.test}}) with
#'constant, linearly modeled, and freely varying scale parameters,
#'respectively.
#'@param exact logical, should an exact permutation test be performed. Only an
#'exact test can be performed for "SO". The default is to use the asymptotic
#'p-values except for "SO".
#'@param R integer, number of permutations for the exact test
#'@param control an optional list of control settings for the stochastic order
#'("SO") test, usually a call to \code{\link{soControl}}.  See there for the
#'names of the settable control values and their effect.
#'@return A list with two components and an optional "boot" attribute that
#'contains the detailed results of the permutation test as an object of class
#'\code{\link[boot]{boot}} if an exact test was performed.
#'@return \item{statistic}{numeric, the value of the test statistic}
#'@return \item{p.val}{numeric, asymptotic one-sided p-value of the test}
#'@author Aniko Szabo, aszabo@@mcw.edu
#'@seealso \code{\link{SO.trend.test}}, \code{\link{RS.trend.test}}, and
#'\code{\link{GEE.trend.test}} for details about the available tests.
#'@keywords htest nonparametric
#'@examples
#'
#'data(shelltox)
#'trend.test(shelltox, test="RS")
#'set.seed(5724)
#'#R=50 is too low to get a good estimate of the p-value
#'trend.test(shelltox, test="RS", R=50, exact=TRUE)
#'

trend.test <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO", 
                       R=100, control=soControl()){ 
   test <- match.arg(test)
   if (!exact && !(test=="SO")){
     res <- switch(test, RS=RS.trend.test(cbdata), 
                         GEE=GEE.trend.test(cbdata,scale.method="fixed"),
                         GEEtrend=GEE.trend.test(cbdata,scale.method="trend"),
                         GEEall=GEE.trend.test(cbdata,scale.method="all"))
   }
   else {
     dat2 <- cbdata[rep(1:nrow(cbdata), cbdata$Freq),]  #each row is one sample
     dat2$Freq <- NULL
     
     boot.LRT.fun <- function(dat, idx){
       dat.new <- cbind(dat[idx, c("ClusterSize","NResp")], Trt=dat$Trt)   #rearrange clusters
       dat.f <- aggregate(dat.new$Trt, 
                  list(Trt=dat.new$Trt, ClusterSize=dat.new$ClusterSize, NResp=dat.new$NResp), length)
       names(dat.f)[names(dat.f)=="x"] <- "Freq"
       dat.f$ClusterSize <- as.numeric(as.character(dat.f$ClusterSize))
       dat.f$NResp <- as.numeric(as.character(dat.f$NResp))
       class(dat.f) <- c("CBData", class(dat.f))
                     
       stat <- switch(test, SO=SO.LRT(dat.f, control=control),
                            RS=RS.trend.test(dat.f)$statistic,
                            GEE=GEE.trend.test(dat.f, scale.method="fixed")$statistic,
                            GEEtrend=GEE.trend.test(cbdata,scale.method="trend")$statistic,
                            GEEall=GEE.trend.test(cbdata,scale.method="all")$statistic)
       stat}        
          
     bootres <- boot(dat2, boot.LRT.fun, R=R, sim="permutation")
     res <- list(statistic=bootres$t0, p.val= mean(bootres$t[,1] >= bootres$t0))
     attr(res, "boot") <- bootres
   }   
   res}

#'Finding the NOSTASOT dose
#'
#'The NOSTASOT dose is the No-Statistical-Significance-Of-Trend dose -- the
#'largest dose at which no trend in the rate of response has been observed. It
#'is often used to determine a safe dosage level for a potentially toxic
#'compound.
#'
#'A series of hypotheses about the presence of an increasing trend overall,
#'with all but the last group, all but the last two groups, etc.  are tested.
#'Since this set of hypotheses forms a closed family, one can test these
#'hypotheses in a step-down manner with the same \code{sig.level} type I error
#'rate at each step and still control the family-wise error rate.
#'
#'The NOSTASOT dose is the largest dose at which the trend is not statistically
#'significant. If the trend test is not significant with all the groups
#'included, the largest dose is the NOSTASOT dose. If the testing sequence goes
#'down all the way to two groups, and a significant trend is still detected,
#'the lowest dose is the NOSTASOT dose. This assumes that the lowest dose is a
#'control group, and this convention might not be meaningful otherwise.
#'
#'@export
#'@param cbdata a \code{\link{CBData}} object
#'@param test character string defining the desired test statistic. See
#'\code{\link{trend.test}} for details.
#'@param exact logical, should an exact permutation test be performed. See
#'\code{\link{trend.test}} for details.
#'@param R integer, number of permutations for the exact test
#'@param sig.level numeric between 0 and 1, significance level of the test
#'@param control an optional list of control settings for the stochastic order
#'("SO") test, usually a call to \code{\link{soControl}}.  See there for the
#'names of the settable control values and their effect.
#'@return a list with two components
#'@return \item{NOSTASOT}{character string identifying the NOSTASOT dose.}
#'@return \item{p}{numeric vector of the p-values of the tests actually performed.}
#'The last element corresponds to all doses included, and will not be missing.
#'p-values for tests that were not actually performed due to the procedure
#'stopping are set to NA.
#'@author Aniko Szabo, aszabo@@mcw.edu
#'@seealso \code{\link{trend.test}} for details about the available trend
#'tests.
#'@references Tukey, J. W.; Ciminera, J. L. & Heyse, J. F. (1985) Testing the
#'statistical certainty of a response to increasing doses of a drug.
#'\emph{Biometrics} 41, 295-301.
#'@keywords htest nonparametric
#'@examples
#'
#'data(shelltox)
#'NOSTASOT(shelltox, test="RS")
#'

NOSTASOT <- function(cbdata, test=c("RS","GEE","GEEtrend","GEEall","SO"), exact=test=="SO",
                     R=100, sig.level=0.05, control=soControl()){
   ntrt <- nlevels(cbdata$Trt)
   control.gr <- levels(cbdata$Trt)[1]
   p.vec <- array(NA, ntrt-1)
   names(p.vec) <- levels(cbdata$Trt)[-1]
   NOSTASOT.found <- FALSE
   curr.gr.idx <- ntrt
   curr.gr <- levels(cbdata$Trt)[ntrt]
   
   while (!NOSTASOT.found & (curr.gr.idx>1)){
     d1 <- cbdata[unclass(cbdata$Trt)<=curr.gr.idx, ]
     d1$Trt <- factor(d1$Trt) #eliminate unused levels
     tr.test <- trend.test(d1, test=test, exact=exact, R=R, control=control)
     p.vec[curr.gr] <- tr.test$p.val
     if (tr.test$p.val < sig.level){ #NOSTASOT not found yet
       curr.gr.idx <- curr.gr.idx - 1
       curr.gr <- levels(cbdata$Trt)[curr.gr.idx]
     }
     else { #NOSTASOT
       NOSTASOT.found <- TRUE
     }
   }
       
   list(NOSTASOT = curr.gr, p=p.vec)    
}       
