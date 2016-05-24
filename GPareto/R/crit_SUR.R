##' Computes the SUR criterion (Expected Excursion Volume Reduction) at point \code{x} for 2 or 3 objectives.
##' To avoid numerical instabilities, the new point is penalized if it is too close to an existing observation.
##' @title Analytical expression of the SUR criterion for two or three objectives.
##' @param x a vector representing the input for which one wishes to calculate the criterion,
##' @param model a list of objects of class \code{\link[DiceKriging]{km}} (one for each objective),
##' @param paretoFront (optional) matrix corresponding to the Pareto front of size \code{[n.pareto x n.obj]}, 
##' @param critcontrol list with two possible options.\cr
##'  
##' A) One can use the four following arguments:
##'        \itemize{
##'        \item \code{integration.points}, matrix of integration points of size \code{[n.integ.pts x d]};
##'        \item \code{integration.weights}, vector of integration weights of length n.integ.pts;
##'        \item \code{mn.X} and \code{sn.X}, matrices of kriging means and sd, each of size \code{[n.obj x n.integ.pts]}; 
##'        \item \code{precalc.data}, list of precalculated data (based on kriging models at integration points) for faster computation. 
##'        }
##' B) Alternatively, one can define arguments passed to \code{\link[GPareto]{integration_design_optim}}:
##'  \code{SURcontrol} (optional), \code{lower}, \code{upper}, \code{min.prob} (optional). This is slower since arguments of
##'  A), used in the function, are then recomputed each time (note that this is not the case when called from \code{\link[GPareto]{GParetoptim}} and \code{\link[GPareto]{crit_optimizer}}).\cr \cr
##' Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} (\code{1e-4)} and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param type "\code{SK}" or "\code{UK}" (default), depending whether uncertainty related to trend estimation has to be taken into account.
##' @seealso \code{\link[GPareto]{crit_EHI}}, \code{\link[GPareto]{crit_SMS}},  \code{\link[GPareto]{crit_EMI}}.
##' @importFrom pbivnorm pbivnorm
##' @importFrom randtoolbox sobol 
##' @return Value of the criterion.
##' @references
##' V. Picheny (2014), Multiobjective optimization using Gaussian process emulators via stepwise uncertainty reduction, 
##' \emph{Statistics and Computing}.
##' @examples
##' #---------------------------------------------------------------------------
##' # crit_SUR surface associated with the "P1" problem at a 15 points design
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' library(DiceDesign)
##' library(KrigInv)
##' 
##' n_var <- 2 
##' n.obj <- 2 
##' f_name <- "P1" 
##' n.grid <- 14
##' test.grid <- expand.grid(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid))
##' n_appr <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
##' response.grid <- t(apply(design.grid, 1, f_name))
##' paretoFront <- t(nondominated_points(t(response.grid)))
##' mf1 <- km(~., design = design.grid, response = response.grid[,1])
##' mf2 <- km(~., design = design.grid, response = response.grid[,2])
##' 
##' model <- list(mf1, mf2)
##' 
##' 
##' integration.param <- integration_design_optim(lower = c(0, 0), upper = c(1, 1), model = model)
##' integration.points <- as.matrix(integration.param$integration.points)
##' integration.weights <- integration.param$integration.weights
##' 
##' precalc.data <- list()
##' mn.X <- sn.X <- matrix(0, nrow = n.obj, ncol = nrow(integration.points))
##' 
##' for (i in 1:n.obj){
##'   p.tst.all <- predict(model[[i]], newdata = integration.points, type = "UK", checkNames = FALSE)
##'   mn.X[i,] <- p.tst.all$mean
##'   sn.X[i,]   <- p.tst.all$sd
##'   precalc.data[[i]] <- precomputeUpdateData(model[[i]], integration.points)
##' }
##' 
##' critcontrol <- list(integration.points = integration.points,
##'                     integration.weights = integration.weights,
##'                     mn.X = mn.X, sn.X = sn.X, precalc.data = precalc.data)
##' ## Alternatively: critcontrol <- list(lower = rep(0, n_var), upper = rep(1,n_var))
##'                 
##' EEV_grid <- apply(test.grid, 1, crit_SUR, model = model, paretoFront = paretoFront,
##'                   critcontrol = critcontrol)
##'                   
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                matrix(pmax(0,EEV_grid), nrow = n.grid), main = "EEV criterion",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             }
##'               )         
##' @export

crit_SUR <- function(x, model, paretoFront = NULL,
                     critcontrol = list(SURcontrol = NULL), type = "UK"){
  X.new <- matrix(x, nrow=1, ncol=model[[1]]@d)
  if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    crit <- -1
  } else {
    

    n.obj <- length(model)
    
    ## Use of SURcontrol for integration_design_optim or directly with integration.points/weights and
    ## mn.X, sn.X
    if(is.null(critcontrol$integration.points)){
      d <- model[[1]]@d
      integration.param   <- integration_design_optim(critcontrol$SURcontrol, d,
                                                      critcontrol$lower, critcontrol$upper, model = model)
      integration.points  <- as.matrix(integration.param$integration.points)
      integration.weights <- integration.param$integration.weights
      
      precalc.data <- vector("list", n.obj)
      mn.X <- sn.X <- matrix(0, n.obj, nrow(integration.points))
      
      for (i in 1:n.obj){
        p.tst <- predict(model[[i]], newdata=integration.points, type=type, checkNames=FALSE)
        mn.X[i,] <- p.tst$mean
        sn.X[i,]   <- p.tst$sd
        if (max(sn.X[i,]) != 0) precalc.data[[i]] <- precomputeUpdateData(model[[i]], integration.points)
      }
      
    }else{
      integration.points  <- critcontrol$integration.points
      integration.weights <- critcontrol$integration.weights
      mn.X <- critcontrol$mn.X 
      sn.X <- critcontrol$sn.X 
      precalc.data <- critcontrol$precalc.data
    }
    
    n.integration.points <- nrow(integration.points)
    if (is.null(integration.weights)) {integration.weights <- rep(1/n.integration.points, n.integration.points)}
    
    if(is.null(paretoFront)){
      observations <- c()    
      for (i in 1:n.obj) observations <- cbind(observations, model[[i]]@y)
      paretoFront <- t(nondominated_points(t(observations)))
    }
    
    #--- Remove integration points with zero variance -----------------------
    I <- which(sn.X[1,] < sqrt(model[[1]]@covariance@sd2)/1e4)
    if (length(I) > 0){
      integration.points <- integration.points[-I,,drop=FALSE]
      mn.X  <- mn.X[,-I, drop=FALSE]
      sn.X  <- sn.X[,-I, drop=FALSE]
      
      for (i in 1:n.obj){
#         if(class(model[[i]]) == "km"){
          precalc.data[[i]]$Kinv.c.olddata <- precalc.data[[i]]$Kinv.c.olddata[,-I,drop=FALSE]
          precalc.data[[i]]$first.member   <- precalc.data[[i]]$first.member[-I]
#         }
      }
      
      if (length(integration.weights) > 1){integration.weights <- integration.weights[-I]}
      
    }
    if (length(integration.weights) < 1) {
      stop("Unable to compute the SUR criterion (crit_SUR): all the integration points have a too small kriging variance")
      pred.xnew <- NULL
      crit <- 0
    } else {
      n.integration.points <- nrow(integration.points)
      integration.weights <- integration.weights / sum(integration.weights)
      
      ##-------- REORDERING --------##
      if (n.obj==2){
        ##-------- 2D: ascending order for obj1 ---------------------##
        if (is.unsorted(paretoFront[,1])) paretoFront <- paretoFront[sort(paretoFront[,1], index.return=TRUE)[[2]],,drop=FALSE]
      } else {
        ##-------- 3D: decreasing order for obj3 --------------------##
        if (is.unsorted(-paretoFront[,3])) paretoFront <- paretoFront[sort(paretoFront[,3], decreasing=TRUE, index.return=TRUE)[[2]],,drop=FALSE]
      }
      n.pareto <- nrow(paretoFront)
      
      ##-------- Precalculations 1 --------##
      mn.xnew <- sn.xnew <- rep(0, n.obj)
      kn <- rho <- eta <- r <- matrix(0, nrow=n.obj, ncol=n.integration.points)
      
      pred.xnew <- vector("list",n.obj)
      for (i in 1:n.obj){
        krig  <- predict(object=model[[i]], newdata=data.frame(x=(X.new)), type=type, se.compute=TRUE,cov.compute=FALSE,checkNames=FALSE)
        mn.xnew[i] <- krig$mean
        sn.xnew[i] <- krig$sd
        pred.xnew[[i]] <- krig
        if (krig$sd != 0) {
          kn[i,] = computeQuickKrigcov2(model[[i]],integration.points=integration.points,X.new=(X.new),
                                        precalc.data=precalc.data[[i]], F.newdata=krig$F.newdata, c.newdata=krig$c)
          rho[i,] <- kn[i,] / (sn.xnew[i]*sn.X[i,])
          denom   <- sqrt( pmax(1e-12, sn.xnew[i]^2 + sn.X[i,]^2 - 2*kn[i,] ))
          eta[i,] <- (mn.xnew[i] - mn.X[i,]) / denom
          r[i,]   <- (kn[i,] - sn.xnew[i]^2) / sn.xnew[i] / denom
        }
      }
      phi.eta <- pnorm( eta )
      
      ##-------- Precalculations for x --------##
      objx.bar <- (paretoFront[,1]-mn.xnew[1]) / sn.xnew[1]
      
      objx.tilde <- (matrix(rep(paretoFront[,1],n.integration.points), ncol = n.integration.points) - t(matrix(rep(mn.X[1,], n.pareto), ncol=n.pareto)) ) / 
        t(matrix(rep(sn.X[1,], n.pareto), ncol=n.pareto))
      
      phi.x.bar   <- pnorm(objx.bar)
      phi.x.tilde <- pnorm(objx.tilde)
      phi.eta.x   <- pnorm(eta[1,])
      
      phi2.x.x <- phi2.x.eta <- matrix(0, (n.pareto), n.integration.points)
      for (i in 1:(n.pareto)){
        phi2.x.x[i,]   <- pbivnorm(rep(objx.bar[i], n.integration.points), objx.tilde[i,],  rho[1,])
        phi2.x.eta[i,] <- pbivnorm(rep(objx.bar[i], n.integration.points), eta[1,], r[1,])
      }
      
      ##-------- Precalculations for y --------##
      if (sn.xnew[2] != 0) {
        objy.bar <- (paretoFront[,2]-mn.xnew[2]) / sn.xnew[2]
        objy.tilde <- (matrix(rep(paretoFront[,2],n.integration.points), ncol = n.integration.points) - t(matrix(rep(mn.X[2,], n.pareto), ncol=n.pareto)) ) / 
          t(matrix(rep(sn.X[2,], n.pareto), ncol=n.pareto))
        
        phi.y.bar   <- pnorm(objy.bar)
        phi.y.tilde <- pnorm(objy.tilde)
        phi.eta.y   <- pnorm(eta[2,])
        
        phi2.y.y <- phi2.y.eta <- matrix(0, (n.pareto), n.integration.points)
        for (i in 1:(n.pareto)) {
          phi2.y.y[i,]   <- pbivnorm(rep(objy.bar[i], n.integration.points), objy.tilde[i,],  rho[2,])
          phi2.y.eta[i,] <- pbivnorm(rep(objy.bar[i], n.integration.points), eta[2,], r[2,])
        }
      } else {
        phi.y.bar   <- as.numeric(paretoFront[,2]>mn.xnew[2])
        phi.y.tilde <- 1*(matrix(rep(paretoFront[,2],n.integration.points), ncol = n.integration.points) > 
                            t(matrix(rep(mn.X[2,], n.pareto), ncol=n.pareto)))
        phi.eta.y   <- as.numeric(mn.xnew[2] > mn.X[2,])
        phi2.y.y <- phi2.y.eta <- matrix(0, (n.pareto), n.integration.points)
        for (i in 1:(n.pareto)) {
          phi2.y.y[i,]   <- phi.y.bar[i]*phi.y.tilde[i,]
          phi2.y.eta[i,] <- phi.y.bar[i]*phi.eta.y
        }
      }
      
      if (n.obj ==3) {
        if (sn.xnew[3] != 0) {
          objz.bar <- (paretoFront[,3]-mn.xnew[3]) / sn.xnew[3]
          objz.tilde <- (matrix(rep(paretoFront[,3],n.integration.points), ncol = n.integration.points) - t(matrix(rep(mn.X[3,], n.pareto), ncol=n.pareto)) ) / 
            t(matrix(rep(sn.X[3,], n.pareto), ncol=n.pareto))
          phi.z.bar   <- pnorm(objz.bar)
          phi.z.tilde <- pnorm(objz.tilde)
          phi.eta.z   <- pnorm(eta[3,])
          
          phi2.z.z <- phi2.z.eta <- matrix(0, (n.pareto), n.integration.points)
          
          for (i in 1:(n.pareto)){
            phi2.z.z[i,]   <- pbivnorm(rep(objz.bar[i], n.integration.points), objz.tilde[i,],  rho[3,])
            phi2.z.eta[i,] <- pbivnorm(rep(objz.bar[i], n.integration.points), eta[3,], r[3,])
          }
        } else {
          phi.z.bar   <- as.numeric(paretoFront[,3]>mn.xnew[3])
          phi.z.tilde <- 1*(matrix(rep(paretoFront[,3],n.integration.points), ncol = n.integration.points) > 
                              t(matrix(rep(mn.X[3,], n.pareto), ncol=n.pareto)))
          phi.eta.z   <- as.numeric(mn.xnew[3] > mn.X[3,])
          
          phi2.z.z <- phi2.z.eta <- matrix(0, (n.pareto), n.integration.points)
          for (i in 1:(n.pareto)){
            phi2.z.z[i,]   <- phi.z.bar[i]*phi.z.tilde[i,]
            phi2.z.eta[i,] <- phi.z.bar[i]*phi.eta.z
          }
        }
      }
      
      if (n.obj==2){
        ######################################################################################################################
        ######## 2D CASE #####################################################################################################
        ######################################################################################################################
        
        res <- EEV.2D.computation(phi.x.bar, phi.x.tilde, phi.eta.x, phi.y.bar, phi.y.tilde, phi.eta.y, 
                                  phi2.x.x, phi2.x.eta, phi2.y.y, phi2.y.eta)
        
        piold <- colSums(res[[2]])
        pinew <- colSums(res[[3]])
        pinew[is.na(pinew)] <- piold[is.na(pinew)]
        
        newpn <- sum(pinew*integration.weights)
        oldpn <- sum(piold*integration.weights)
        
        crit <- oldpn - newpn
      } else
      {
        ######################################################################################################################
        ######## 3D CASE #####################################################################################################
        ######################################################################################################################
        
        #*** Current excursion volume **********************************************
        pn <- pold <- rep(0, n.integration.points)
        
        #*** START new excursion volume ********************************************
        for (i in 1:(n.pareto)){
          
          nondominated.sub <- which(!is_dominated((nondominated_points(t(paretoFront[(i:n.pareto),1:2, drop=FALSE])))))
          #         nondominated.sub <- paretocheck(paretoFront[(i:n.pareto),1:2, drop=FALSE])                         # Subset of the Pareto front
          nondominated.sub <- i-1+nondominated.sub[sort(paretoFront[(i:n.pareto)[nondominated.sub],1],index.return=TRUE)[[2]]] # Reorder by increasing first objective
          n.pareto.sub     <- length(nondominated.sub)
          
          res <- EEV.2D.computation(phi.x.bar[nondominated.sub], phi.x.tilde[nondominated.sub,,drop=FALSE], phi.eta.x, 
                                    phi.y.bar[nondominated.sub], phi.y.tilde[nondominated.sub,,drop=FALSE], phi.eta.y, 
                                    phi2.x.x[nondominated.sub,,drop=FALSE], phi2.x.eta[nondominated.sub,,drop=FALSE], 
                                    phi2.y.y[nondominated.sub,,drop=FALSE], phi2.y.eta[nondominated.sub,,drop=FALSE])
          
          pijold <- colSums(res[[2]])
          pijnew <- colSums(res[[3]])
          
          if (i==1) {
            p.joint.1 <- 1 - phi2.z.z[i,] - phi.eta.z + phi2.z.eta[i,]
            p.joint.2 <- phi2.z.z[i,] - phi.z.tilde[i,] + phi.eta.z - phi2.z.eta[i,]
          } else {
            p.joint.1 <- phi2.z.z[i-1,] - phi2.z.z[i,] - phi2.z.eta[i-1,] + phi2.z.eta[i,]
            p.joint.2 <- phi.z.tilde[i-1,] - phi.z.tilde[i,] - phi2.z.z[i-1,] + phi2.z.z[i,] + phi2.z.eta[i-1,] - phi2.z.eta[i,]          
          }
          pold <- pold + pijold*(p.joint.1+p.joint.2)
          pn   <- pn + pijnew*(p.joint.1) + pijold*(p.joint.2)
        }
        # Add last slice
        pold <- pold + phi.z.tilde[n.pareto,]
        pa <- phi.z.tilde[n.pareto,] - phi2.z.z[n.pareto,] + phi2.z.eta[n.pareto,]
        px <- phi.eta.x
        py <- phi.eta.y
        pb <- ( phi2.z.z[n.pareto,] - phi2.z.eta[n.pareto,]) * (px + py - px*py )
        pn <- pn + (pa + pb)
        pn[is.na(pn)] <- pold[is.na(pn)]
        newpn <- sum(pn*integration.weights)
        oldpn <- sum(pold*integration.weights)
        #*** STOP new excursion volume *********************************************
        crit <- oldpn - newpn
      }
    }
    if (crit<1e-32) crit <- prob.of.non.domination(paretoFront=paretoFront, model, X.new, predictions=pred.xnew)-1

    return(crit)
  }
  
  #   return(list(crit=(1-crit*n.integration.points/(n.integration.points+1))) )
}