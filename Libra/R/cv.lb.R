#' CV for lb
#' 
#' Cross-validation method to tuning the parameter t for lb.
#' 
#' K-fold cross-validation method is used to tuning the parameter t for ISS.
#' Mean square error is used for linear model. Miss-classification error
#' is used for binomial and multinomial model.
#' 
#' @param X An n-by-p matrix of predictors
#' @param y Response Variable
#' @param kappa The damping factor of the Linearized Bregman Algorithm that is
#'  defined in the reference paper. See details. 
#' @param alpha Parameter in Linearized Bregman algorithm which controls the 
#' step-length of the discretized solver for the Bregman Inverse Scale Space. 
#' See details. 
#' @param K Folds number for CV. Default is 5.
#' @param tlist Parameters t along the path.
#' @param nt Number of t. Used only if tlist is missing. Default is 100.
#' @param trate tmax/tmin. Used only if tlist is missing. Default is 100.
#' @param family Response type
#' @param group Whether to use a group penalty, Default is FALSE.
#' @param intercept If TRUE, an intercept is included in the model (and not 
#' penalized), otherwise no intercept is included. Default is TRUE.
#' @param normalize if TRUE, each variable is scaled to have L2 norm 
#' square-root n. Default is TRUE.
#' @param plot.it Plot it? Default is TRUE
#' @param se Include standard error bands? Default is TRUE
#' @param \dots Additonal arguments passing to lb
#' @return A list is returned. The list contains a vector of parameter t, 
#' crossvalidation error cv.error, and the estimated standard deviation for it cv.sd
#' @references Ohser, Ruan, Xiong, Yao and Yin, Sparse Recovery via Differential
#'  Inclusions, \url{http://arxiv.org/abs/1406.7728}
#' @author Feng Ruan, Jiechao Xiong and Yuan Yao
#' @keywords Cross-validation
#' @examples
#' #Examples in the reference paper
#' library(MASS)
#' n = 200;p = 100;k = 30;sigma = 1
#' Sigma = 1/(3*p)*matrix(rep(1,p^2),p,p)
#' diag(Sigma) = 1
#' A = mvrnorm(n, rep(0, p), Sigma)
#' u_ref = rep(0,p)
#' supp_ref = 1:k
#' u_ref[supp_ref] = rnorm(k)
#' u_ref[supp_ref] = u_ref[supp_ref]+sign(u_ref[supp_ref])
#' b = as.vector(A%*%u_ref + sigma*rnorm(n))
#' cv.lb(A,b,10,1/20,intercept = FALSE,normalize = FALSE)
#' 
#' #Simulated data, binomial case
#' X <- matrix(rnorm(500*100), nrow=500, ncol=100)
#' alpha <- c(rep(1,30), rep(0,70))
#' y <- 2*as.numeric(runif(500)<1/(1+exp(-X %*% alpha)))-1
#' cv.lb(X,y,kappa=5,alpha=1,family="binomial",
#'              intercept=FALSE,normalize = FALSE)
#' 

cv.lb <-function(X, y,kappa,alpha,K = 5, tlist,nt = 100,trate = 100, family = c("gaussian", "binomial", "multinomial"),
                 group = FALSE,intercept = TRUE,normalize = TRUE,plot.it = TRUE, se = TRUE,...)
{
  family <- match.arg(family)
  if (!is.matrix(X)) stop("X must be a matrix!")
  if (family!="multinomial"){
    if (!is.vector(y)) stop("y must be a vector unless in multinomial model!")
    if (nrow(X) != length(y)) stop("Number of rows of X must equal the length of y!")
    if (family=="binomial" & any(abs(y)!=1)) stop("y must be in {1,-1}")
  }
  if (family=="multinomial"){
    if (is.vector(y)){
      if(nrow(X) != length(y)) stop("Number of rows of X must equal the length of y!")
      y_unique <- unique(y)
      y = sapply(1:length(y_unique),function(x) as.numeric(y==y_unique[x]))
    }
    else if (is.matrix(y)){
      if(nrow(X) != nrow(y)) stop("Number of rows of X and y must equal!")
      if (any((y!=1)&(y!=0)) || any(rowSums(y)!=1)) stop("y should be indicator matrix!")
    }
    else
      stop("y must be a vector or matrix!")
  }
  if (group){
    if (missing(index)){
      if (family=="multinomial"){
        index=NA
      }else{
        group=FALSE
        print("Index is missing, using group=FALSE instead!")
      }
    }
    if (!is.vector(index)) stop("Index must be a vector!")
    if ((length(index) + intercept) != ncol(X))
      stop("Length of index must be the same as the number of columns of X minus the intercept!")
  }
  
  n <- dim(X)[1]
  folds <- split(sample(1:n), rep(1:K, length = n))
  if(missing(tlist)){
    obj <- lb(X, y, kappa, alpha, nt=nt,trate=trate, family = family, group=group,
              intercept = intercept, normalize=normalize,...)
    tlist <- obj$t
  }
  
  residmat <- matrix(0, length(tlist), K)
  for(i in seq(K)) {
    omit <- folds[[i]]
    if (family=="multinomial")
        fit <- lb(X[-omit,  ,drop=FALSE], y[-omit,  ,drop=FALSE], kappa, alpha, tlist,nt=nt,trate=trate,family = family, group=group,
                intercept = intercept, normalize=normalize,...)
    else
        fit <- lb(X[-omit,  ,drop=FALSE], y[-omit,drop=FALSE], kappa, alpha,tlist,nt=nt,trate=trate, family = family, group=group,
                intercept = intercept, normalize=normalize,...)
    fit <- predict(fit, X[omit,  ,drop=FALSE], tlist)$fit
    
    if (family=="multinomial")
        residmat[,i] <- sapply(1:length(t), function(x) sum((y[omit,,drop=FALSE]-(fit[[x]]==apply(fit[[x]],1,max)))^2)/2/length(omit))
    else if(length(omit)==1)fit<-matrix(fit,nrow=1)
    if (family=="binomial")
        residmat[, i] <- apply((y[omit] - 2*(fit>0.5)+1)^2/4, 2, mean)
    if (family=="gaussian")
        residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
  }
  cv.error <- apply(residmat, 1, mean)
  cv.sd <- sqrt(apply(residmat, 1, var)/K)
  object<-list(t = tlist, cv.erroe = cv.error, cv.sd = cv.sd)
  if (plot.it){
    plot(tlist, cv.error,pch="*",type = "b", ylim = range(cv.error, cv.error + cv.sd, 
                                                      cv.error - cv.sd), xlab=("t"),ylab="Cross-Validated MSE")
    if(se)
      segments(tlist, cv.error + cv.sd, tlist, cv.error - cv.sd, col="red")
  }
  object
}

