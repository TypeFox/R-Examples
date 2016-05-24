#' @title Sub-window permutation analysis coupled with PLS (SwPA-PLS)
#'
#' @description SwPA-PLS provides the influence of each variable without considering the 
#' influence of the rest of the variables through sub-sampling of samples and variables.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param N number of Monte Carlo simulations (default = 3).
#' @param ratio the proportion of the samples to use for calibration (default = 0.8).
#' @param Qv integer number of variables to be sampled in each iteration (default = 10).
#' @param SPA.threshold thresholding to remove non-important variables (default = 0.05).
#'  
#' @return Returns a vector of variable numbers corresponding to the model 
#' having lowest prediction error.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references H. Li, M. Zeng, B. Tan, Y. Liang, Q. Xu, D. Cao, Recipe for revealing 
#' informative metabolites based on model population analysis, Metabolomics 6 (2010) 353-361.
#' http://code.google.com/p/spa2010/downloads/list.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' with( gasoline, spa_pls(octane, NIR) )
#'
#' @export
spa_pls<- function(y, X, ncomp=10, N=3, ratio=0.8, Qv=10, SPA.threshold=0.05){
  
  # Subwindow Permutation Analysis for variable assessment.
  #	   Input:  X: m x p  (Sample matrix)
  #            y: m x 1  (measured property)
  #      ncomp: The allowed maximum number of PLS components for cross-validation
  #            N: The number of Monte Carlo Simulation.
  #        ratio: The ratio of calibration samples to the total samples.
  #           Qv: The number of variables to be sampled in each MCS.
  
  #  Output: Structural data: F with items:
  #        model: a matrix of size N x p with element 0 or 1(means the variable is selected).
  #          nLV: The optimal number of PLS components for each submodel.
  #       error0: The normal prediction errors of the N submodels.
  #       error1: The permutation prediction errors of the N submodels
  #     interfer: a vector of size p. '1' indicates the variable is interferring
  #            p: the p-value of each variable resulting from SPA. 
  
  modeltype <- "prediction"
  if (is.factor(y)) {
    modeltype <- "classification"
    tb<-as.numeric(names(table(y)))
  }
  
  Mx <- dim(X)[1]
  Nx <- dim(X)[2] 
  Qs <- floor(Mx*ratio)
  error0   <- matrix(0,  N,1)
  error1   <- matrix(NA, N,Nx)
  interfer <- matrix(0,  1,Nx)
  nLV <- matrix(0, N,1)
  ntest <- Mx-Qs
  
  for (i in 1:N){
    ns    <- sample(Mx)
    calk  <- ns[1:Qs]
    testk <- ns[(Qs+1):Mx]
    nv <- sample(Nx)[1:Qv]
    variableIndex <- matrix(0,1,Nx)
    variableIndex[nv] <- 1   
    
    Xcal  <- X[calk,nv];  ycal  <- y[calk]
    Xtest <- X[testk,nv]; ytest <- y[testk]    
    
    pls.object <- plsr(ycal ~ Xcal, ncomp=min(ncomp, (ncol(Xcal)-1)), validation = "LOO")
    Press    <- pls.object$valid$PRESS[1,]
    opt.comp <- which.min(Press)
    yy <- c(ycal, ytest); XX <- rbind( Xcal, Xtest)
    mydata  <- data.frame( yy=yy ,XX=I(XX) , train=c(rep(TRUE, length(ycal)), rep(FALSE, length(ytest))))
    pls.fit <- plsr(yy ~ XX, ncomp=opt.comp, data=mydata[mydata$train,])
    
    # Make predictions using the established PLS model 
    if (modeltype == "prediction"){
      pred.pls  <- predict(pls.fit, ncomp = opt.comp, newdata=mydata[!mydata$train,])
      error0[i] <- sqrt(sum((ytest-pred.pls[,,])^2))
    } else if (modeltype == "classification"){
      may.lda   <- lda(y ~ X.score, data = data.frame(y = ycal, X.score = I(pls.fit$scores )))
      score.val <- (Xtest - rep(colMeans(Xcal),each = nrow(Xtest))) %*% pls.fit$projection 
      error0[i] <- 100-sum(ytest == predict(may.lda, dimen = opt.comp, newdata = data.frame(X.score = I(score.val)))$class)/nrow(Xtest)*100
    }
    error_temp <- matrix(NA,1,Nx)
    
    # Make predictions on the permutated sub-datset.
    for (j in 1:Qv){
      rn     <- sample(ntest)
      Xtestr <- Xtest
      vi     <- Xtest[,j]
      Xtestr[,j] <- vi[rn]
      newdata    <- data.frame( yy=ytest ,XX=I(Xtestr))
      if (modeltype == "prediction"){
        predr.pls <- predict(pls.fit, ncomp = 1:opt.comp, newdata=newdata)
        error_temp[nv[j]] <- sqrt(sum((ytest-predr.pls[,,])^2))
      } else if (modeltype == "classification"){
        may.lda   <- lda(y ~ X.score, data = data.frame( y = ycal, X.score = I(pls.fit$scores )))
        score.val <- (Xtestr - rep(colMeans(Xcal),each = nrow(Xtestr))) %*% pls.fit$projection 
        error_temp[nv[j]] <- 100-sum(ytest == predict(may.lda, dimen = opt.comp, newdata = data.frame(X.score = I(score.val)))$class)/nrow(Xtestr)*100
      }
    }
    error1[i,] <- error_temp
  }
  #p value computing
  p     <- matrix(0, 1,Nx)
  DMEAN <- matrix(0, Nx,1)
  DSD   <- matrix(0, Nx,1)
  for (ii in 1:Nx){
    k <- which(!is.na(error1[,ii]))
    if(length(k)!=0){
      errori <- error1[,ii]
      errorn <- error0[k]
      errorp <- errori[k]
      MEANn  <- mean(errorn)
      MEANp  <- mean(errorp)
      SDn <- sd(errorn)
      SDp <- sd(errorp)
      DMEAN[ii] <- MEANp - MEANn
      DSD[ii]   <- SDp - SDn
      ranksum.test <- wilcox.test(errorn,errorp, paired = TRUE, alternative = "greater")
      p[ii] <- ranksum.test$p.value
    }
  }
  
  spa.selection <- which(p < SPA.threshold)
  if(length(spa.selection) <= (ncomp +1)) {
    spa.selection <-sort(p,decreasing = T, index.return = T)$ix [1:ncomp]
  }
  return(list(spa.selection=spa.selection))
}
