##' Calculates GoF measures for Cox's propoportional hazard model for right
##' censored survival times
##' 
##' Calculates score processes and KS and Cvm tests for proportionaly of hazards
##' via simulation (Martinussen and Scheike, 2006).
##' 
##' 
##' @param model Model object (\code{lm} or \code{glm})
##' @param variable List of variable to order the residuals after
##' @param R Number of samples used in simulation
##' @param type Type of GoF-procedure
##' @param plots Number of realizations to save for use in the plot-routine
##' @param seed Random seed
##' @param ... additional arguments
##' @return Returns an object of class 'cumres'.
##' @author Klaus K. Holst and Thomas Scheike
##' @method cumres coxph
##' @export
##' @seealso \code{\link[gof]{cumres.glm}}, \code{\link[survival]{coxph}}, and
##' \code{\link[timereg]{cox.aalen}} in the \code{timereg} package for similar
##' GoF-methods for survival-data.
##' @references Lin, D. Y. and Wei, L. J. and Ying, Z. (1993) \emph{Checking the
##' Cox model with cumulative sums of martingale-based residuals} Biometrika,
##' Volume 80, No 3, p. 557-572.
##' 
##' Martinussen, Torben and Scheike, Thomas H.  \emph{Dynamic regression models
##' for survival data} (2006), Springer, New York.
##' @keywords models regression
##' @examples
##' 
##' library(survival)
##' 
##' simcox <- function(n=100, seed=1) {
##'   if (!is.null(seed))
##'     set.seed(seed)
##'   require(survival)
##'   time<-rexp(n); cen<-2*rexp(n);
##'   status<-(time<cen);
##'   time[status==0]<-cen[status==0];
##'   X<-matrix(rnorm(2*n),n,2)
##'   return(data.frame(time=time, status=status, X))
##' }
##' n <- 100; d <- simcox(n); m1 <- coxph(Surv(time,status)~ X1 + X2, data=d)
##' cumres(m1)
##' 
##' \dontrun{
##' ## PBC example
##' data(pbc)
##' fit.cox <- coxph(Surv(time,status==2) ~ age + edema + bili + protime, data=pbc)
##' system.time(pbc.gof <- cumres(fit.cox,R=2000))
##' par(mfrow=c(2,2))
##' plot(pbc.gof, ci=TRUE, legend=NULL)
##' }
`cumres.coxph` <- function(model,
         variable=c(colnames(model.matrix(model))),
         type=c("score","residual"),
         R=1000, plots=min(R,50), seed=round(runif(1,1,1e9)), ...) {

  require(survival)
  mt <- model.frame(model)
  Y <- model.extract(mt, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  if (attr(Y, "type") == "right") {
    time <- Y[, "time"]; 
    status <- Y[, "status"]
  } else stop("Expected right-censored data.");
##  X <- na.omit(model.matrix(model)[,-1,drop=FALSE]) ## Discard intercept
  X <- na.omit(model.matrix(model))

  ot <- order(time,status==0); # order in time, status=1 first for ties
  time <- time[ot]; status <- status[ot]
  X <- X[ot,,drop=FALSE]
  n <- length(time)
  nd <- sum(status)
  p <- ncol(X)
  index.dtimes <- (1:n)[status==1]
  dtimes <- time[index.dtimes]

  Iinv <- model$naive.var
  beta.iid <- matrix(residuals(model,type="dfbeta"),ncol=p)[ot,,drop=FALSE]
  ##  Mres <- Mt <- residuals(model, type="martingale")[ot]
  ##  cox.schoen <- residuals(model,type="schoen")[,ot]
  beta <- coef(model)
  if(any(is.na(beta))) stop("Over-parametrized model")

  ##  if (!is.na(match(response, variable))) variable[match(response, variable)] <- "predicted"
  if (is.numeric(variable))
    variable <- na.omit(colnames(X)[variable])
  variable <- unique(variable)
  UsedData <- X[,na.omit(match(variable, colnames(X))),drop=FALSE]
  
  myvars <- colnames(UsedData)
  ## Martingale residuals only cumulated after variables with more than two levels
  ##  myvars <- colnames(UsedData)[apply(UsedData,2,function(x) length(unique(x))>2)] ## Only consider variables with more than two levels
  ##  if ("predicted"%in%variable) myvars <- c("predicted",myvars)
  myvars.idx <- which(colnames(X) %in% myvars)

  hatW.MC <- function(x) {
    output <- .C("coxscoreW",
                 R=as.integer(R),
                 n=as.integer(n),
                 nd=as.integer(nd),
                 p=as.integer(p),
                 beta_data=as.double(beta),
                 time_data=as.double(time),
                 index_dtimes_data=as.integer(index.dtimes-1),
                 X_data=as.double(X), # nxp
                 beta_iid_data=as.double(beta.iid), # nxp
                 Mt_data=as.double(as.numeric(n)),
                 paridx=as.integer(x-1), nparidx=as.integer(1),
                 Type=as.integer(1), # 1=Score process
                 seed=round(runif(1,1,1e9)),
                 plotnum=as.integer(plots),
                 
                 KS=as.double(0),
                 CvM=as.double(0),
                 Wsd=as.double(numeric(nd)),
                 cvalues=as.double(numeric(R)),
                 Ws=as.double(numeric(nd*plots)),
                 W=as.double(numeric(nd)),
                 WWW=as.double(0), ## Only for debugging
                 pkg="gof"
                 )
    return(output)
  }

  UsedVars <- W <- Wsd <- What <- KS <- CvM <- allcvalues <- x <- mytype <- c()
  ### SCORE:  
  for (i in 1:length(myvars)) {
    UsedVars <- c(UsedVars, myvars[i])
    onesim <- hatW.MC(myvars.idx[i])
    ##onesim$WWW <- matrix(onesim$WWW,ncol=nd,nrow=n, byrow=TRU)E/    onesim$WWW <- matrix(onesim$WWW,ncol=n,nrow=nd,byrow=TRUE)
    ##    x <- cbind(x, sort(X[,myvars[i]])); colnames(x)[ncol(x)] <- myvars[i]
    x <- cbind(x, dtimes)
    W <- cbind(W, onesim$W)
    Wsd <- cbind(Wsd, onesim$Wsd)
    What <- c(What, list(matrix(onesim$Ws,ncol=plots)));
    KS <- c(KS, onesim$KS);  CvM <- c(CvM, onesim$CvM)
    allcvalues <- cbind(allcvalues, onesim$cvalues)
    mytype <- c(mytype, "score")
  }

  
      
  res <- list(W=W, What=What,
              x=x, 
              KS=KS, CvM=CvM,
              cvalues=allcvalues, variable=UsedVars,
              R=R, sd=Wsd, type=mytype, model="coxph")
  class(res) <- "cumres"
  res
}

