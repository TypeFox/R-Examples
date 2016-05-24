#' @title Uninformative variable elimination in PLS (UVE-PLS)
#'
#' @description Artificial noise variables are added to the predictor set before the PLSR 
#' model is fitted. All the original variables having lower "importance" than the artificial 
#' noise variables are eliminated before the procedure is repeated until a stop criterion is 
#' reached.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param N number of samples Mone Carlo simulations (default = 3).
#' @param ratio the proportion of the samples to use for calibration (default = 0.75).
#' @param MCUVE.threshold thresholding separate signal from noise (default = NA creates 
#' automatic threshold from data).
#'  
#' @return Returns a vector of variable numbers corresponding to the model 
#' having lowest prediction error.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references V. Centner, D. Massart, O. de Noord, S. de Jong, B. Vandeginste, C. Sterna, 
#' Elimination of uninformative variables for multivariate calibration, Analytical Chemistry 
#' 68 (1996) 3851-3858.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' with( gasoline, mcuve_pls(octane, NIR) )
#'
#' @export
mcuve_pls<- function(y,X,ncomp=10, N=3,ratio=0.75, MCUVE.threshold=NA){
  # 	 uninformative variable elimination(UVE)-PLS.
  #	 Input:    y: m x 1  (measured property)
  #            X: m x n  (Sample matrix)
  #        ncomp: The max PC for cross-validation
  #            N: The number of Monte Carlo Simulation.
  #        ratio: The ratio of calibration samples to the total samples.
  
  Mx <- nrow(X)
  Nx <- ncol(X)
  W  <- matrix(runif(Nx*Mx,0,1), Mx, Nx)
  Z  <- cbind(X, W)
  K  <- floor(Mx*ratio) 
  C  <- matrix(0, N, Nx*2)  
  ncomp <- min(c(dim(X), ncomp))
  for( i in 1:N){
    temp <- sample(Mx)
    calk <- temp[1:K]      
    Zcal <- Z[calk, ]; ycal <- y[calk]  
    pls.object <- plsr(ycal ~ Zcal,  ncomp=min(ncomp, (ncol(Zcal)-1)), validation = "LOO")
    Press      <- pls.object$valid$PRESS[1,]
    opt.comp   <- which.min(Press)
    C[i,] <- pls.object$coefficients[,1,opt.comp]
  }
  U  <- apply(C, 2, mean)  
  SD <- apply(C, 2, sd)  
  RI <- U/SD	
  if(is.na( MCUVE.threshold)) {
    MCUVE.threshold <- max(abs(RI[-(1:Nx)])) # max.RI
  } 
  UVE.selection <- which(abs(RI[(1:Nx)]) > MCUVE.threshold)
  if(length(UVE.selection)<= (ncomp +1)) {
    UVE.selection <- sort(abs(RI[(1:Nx)]), decreasing=TRUE, index.return = T)$ix [1:ncomp]
  }
  return(list(mcuve.selection=UVE.selection))
}
