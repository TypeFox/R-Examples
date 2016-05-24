##' Display results of multi-objective optimization returned by either \code{\link[GPareto]{GParetoptim}} or \code{\link[GPareto]{easyGParetoptim}},
##' possibly completed with various post-processings of uncertainty quantification.  
##' @title Plot multi-objective optimization results and post-processing
##' @param res list returned by \code{\link[GPareto]{GParetoptim}} or \code{\link[GPareto]{easyGParetoptim}},
##' @param add logical; if \code{TRUE} adds the first graphical output to an already existing plot; 
##' if \code{FALSE}, (default) starts a new plot,
##' @param UQ_PF logical; for 2 objectives, if \code{TRUE} perform a quantification of uncertainty on
##'  the Pareto front to display the symmetric deviation function with \code{\link[GPareto]{plotSymDevFun}} (cannot be added to existing graph),
##' @param UQ_PS logical; if \code{TRUE} call \code{\link[GPareto]{plot_uncertainty}} representing the probability of non-domination in the variable space,
##' @param UQ_dens logical; for 2D problems, if \code{TRUE} call \code{\link[GPareto]{ParetoSetDensity}} to estimate and display the density of Pareto optimal points in the variable space,
##' @param lower optional vector of lower bounds for the variables. 
##' Necessary if \code{UQ_PF} and/or \code{UQ_PS} are \code{TRUE} (if not provided, variables are supposed to vary between 0 and 1),
##' @param upper optional vector of upper bounds for the variables. 
##' Necessary if \code{UQ_PF} and/or \code{UQ_PS} are \code{TRUE} (if not provided, variables are supposed to vary between 0 and 1),
##' @param control optional list, see details.
##' @details 
##' By default, \code{plotGPareto} displays the Pareto front delimiting the non-dominated area with 2 objectives,
##'  by a perspective view with 3 objectives and using parallel coordinates with more objectives.\cr 
##' 
##' Setting one or several of UQ_PF, UQ_PS and UQ_dens allows to run and display post-processing tools that assess
##' the precision and confidence of the optimization run, either in the objective (\code{UQ_PF}) or the variable spaces 
##' (\code{UQ_PS}, \code{UQ_dens}). Note that these options are computationally intensive.
##' 
##' Various parameters can be used for the display of results and/or passed to subsequent function:
##' \itemize{
##'   \item \code{col}, \code{pch} correspond the color and plotting character for observations,
##'   \item \code{PF.line.col}, \code{PF.pch}, \code{PF.points.col} define the color of the line denoting the current Pareto front,
##'   the plotting character and color of non-dominated observations, respectively,
##'   \item \code{nsim}, \code{npsim} and \code{gridtype} define the number of conditional simulations performed with \code{\link[DiceKriging]{simulate}}
##'    along with the number of simulation points (in case \code{UQ_PF} and/or \code{UQ_dens} are \code{TRUE}),
##'   \item \code{gridtype} to define how simulation points are selected; 
##'   alternatives are '\code{runif}' (default) for uniformly sampled points,
##'   '\code{LHS}' for a Latin Hypercube design using \code{\link[DiceDesign]{lhsDesign}} and 
##'   '\code{grid2d}' for a two dimensional grid,
##'   \item \code{f1lim}, \code{f2lim} can be passed to \code{\link[GPareto]{CPF}},
##'   \item \code{resolution}, \code{option}, \code{nintegpoints} are to be passed to \code{\link[GPareto]{plot_uncertainty}}
##'   \item \code{displaytype} type of display for \code{UQ_dens}, see \code{\link[ks]{plot.kde}},
##'   \item \code{printVD} logical, if \code{TRUE} and \code{UQ_PF} is \code{TRUE} as well, print the value of the Vorob'ev deviation,
##'   \item \code{meshsize3d} mesh size of the perspective view for 3-objective problems,
##'   \item \code{theta}, \code{phi} angles for perspective view of 3-objective problems.
##' }
##' @export
##' @importFrom graphics abline persp 
##' @references
##' M. Binois, D. Ginsbourger and O. Roustant (2015), Quantifying Uncertainty on Pareto Fronts with Gaussian process conditional simulations, 
##' \emph{European Journal of Operational Research}, 243(2), 386-394. \cr \cr
##' A. Inselberg (2009), \emph{Parallel coordinates}, Springer.
##' @examples 
##' \dontrun{
##' #---------------------------------------------------------------------------
##' # 2D objective function
##' #---------------------------------------------------------------------------
##' set.seed(25468)
##' n_var <- 2 
##' fname <- P1
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' res <- easyGParetoptim(fn=fname, lower=lower, upper=upper, budget=15, 
##' control=list(method="EHI", inneroptim="pso", maxit=20))
##' 
##' ## Pareto front only
##' plotGPareto(res)
##' 
##' ## With post-processing
##' plotGPareto(res, UQ_PF = TRUE, UQ_PS = TRUE, UQ_dens = TRUE)
##' 
##' }
plotGPareto <- function(res, add = FALSE, UQ_PF = FALSE, UQ_PS = FALSE, UQ_dens = FALSE,
                        lower = NULL, upper = NULL,
                        control = list(pch = 20, col = "red", PF.line.col = "cyan", PF.pch = 17,
                                       PF.points.col = "blue", VE.line.col = "cyan",
                                       nsim = 100, npsim = 1500, gridtype = "runif",
                                       displaytype = "persp", printVD = TRUE,
                                       meshsize3d = 50, theta = -25, phi = 10)){
  # Check of arguments
  if(is.null(control$pch)) control$pch <- 20
  if(is.null(control$col)) control$col <- "red"
  if(is.null(control$PF.line.col)) control$PF.line.col <- "cyan"
  if(is.null(control$PF.pch)) control$PF.pch <- 17
  if(is.null(control$VE.line.col)) control$VE.line.col <- "cyan"
  if(is.null(control$PF.points.col)) control$PF.points.col <- "blue"
  if(is.null(control$nsim)) control$nsim <- 100
  if(is.null(control$npsim)) control$npsim <- 1000
  if(is.null(control$gridtype)) control$gridtype <- "runif"
  if(is.null(control$displaytype)) control$displaytype <- "persp"
  if(is.null(control$printVD)) control$printVD <- TRUE
  if(is.null(control$meshsize3d)) control$meshsize3d <- 50
  if(is.null(control$theta)) control$theta <- -25
  if(is.null(control$phi)) control$phi <- 10
  
  if(is.null(lower)) lower <- rep(0, ncol(res$par))
  if(is.null(upper)) upper <- rep(1, ncol(res$par))
  
  n.obj <- ncol(res$value)
  
  if(!UQ_PF){
    if(n.obj == 2){
      if(!add){
        ## easyGParetoptim
        if(!is.null(res$history)){
          plot(res$history$y, pch = control$pch, col = control$col, xlab = expression(f[1]), ylab = expression(f[2])) 
        }else{## GParetoptim
          plot(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y,
               pch = control$pch, col = control$col, xlab = expression(f[1]), ylab = expression(f[2]) )
        }
      }else{
        if(!is.null(res$history)){
          points(res$history$y, pch = control$pch, col = control$col)
        }else{
          points(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y,
                 pch = control$pch, col = control$col)
        }
      }
      
      if(!is.null(res$history)){
        plotParetoEmp(res$value, col = control$PF.line.col)
        points(res$value, col = control$PF.points.col, pch = control$PF.pch)
      }else{
        plotParetoEmp(t(nondominated_points(t(cbind(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y)))),
                      col = control$PF.line.col)
        points(t(nondominated_points(t(cbind(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y)))),
               col = control$PF.points.col, pch = control$PF.pch)
      }
    }else{
      if(n.obj == 3){
        if(!is.null(res$history)){
          ally <- res$value
        }else{
          ally <- cbind(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y, res$lastmodel[[3]]@y)
        }
        
        ally <- apply(ally, c(1,2), round, digits = 4)
        ally <- ally[!duplicated(ally),]
        
        xax <- sort(ally[,1])
        xax <- xax[!duplicated(xax)]
        
        yax <- sort(ally[,2])
        yax <- yax[!duplicated(yax)]
        
        zlevels <- sort(ally[,3])
        zlevels <- zlevels[!duplicated(zlevels)]
        
        # add levels if there are not as many as meshsize3d^3
        if(length(xax) * length(yax) * length(zlevels) < control$meshsize3d^3){
          while(length(xax) < control$meshsize3d){
            tmp <- order(xax[-1] - xax[-length(xax)], decreasing = TRUE)[1]
            xax <- sort(c(xax, (xax[tmp + 1] + xax[tmp])/2))
          }
          while(length(yax) < control$meshsize3d){
            tmp <- order(yax[-1] - yax[-length(yax)], decreasing = TRUE)[1]
            yax <- sort(c(yax, (yax[tmp + 1] + yax[tmp])/2))
          }
          while(length(zlevels) < control$meshsize3d){
            tmp <- order(zlevels[-1] - zlevels[-length(zlevels)], decreasing = TRUE)[1]
            zlevels <- sort(c(zlevels, (zlevels[tmp + 1] + zlevels[tmp])/2))
          }
        }
        
        # rawGrid <- as.matrix(expand.grid(ally[,1], ally[,2], ally[,3]))
        rawGrid <- as.matrix(expand.grid(xax, yax, zlevels))
        dom_elements <- apply(rawGrid, 1, function(x){is_dominated(cbind(as.vector(x), t(ally + 1e-6)))[1]})
        rawGrid <- rawGrid[!dom_elements,]
        dom_elements <- apply(rawGrid, 1, function(x){is_dominated(cbind(as.vector(x), t(ally - 1e-6)))[1]})
        rawGrid <- rawGrid[dom_elements,]
        

        
        xygrid <- as.matrix(expand.grid(xax, yax))
        z <- apply(xygrid, 1, function(x){
                                          tmp <- rawGrid[which(rawGrid[,1] == x[1] & rawGrid[,2] == x[2]),3]
                                          if(length(tmp) == 0) return(NA)
                                          return(min(tmp))
                                          }
        )
        
        persp(x = xax, y = yax, z = matrix(z, nrow = length(xax)), theta = control$theta, phi = control$phi, scale = TRUE,
              ticktype = "detailed", xlab = "f1", ylab = "f2", zlab = "f3")

        
      }else{
        if(!is.null(res$history)){
          parplotPF_nd(res$history$y, add = add)
        }else{
          ally <- NULL
          for (i in 1:length(res$lastmodel)) ally <- cbind(ally, res$lastmodel[[i]]@y)
          parplotPF_nd(ally, add = add)
        }
      }
    }
  }else{
    
    if(n.obj > 2){
      cat("UQ_PF is only implemented for bi-objective problems. \n")
      return(0)
    }
    
    if(control$gridtype == "runif"){
      simPoints <- matrix(runif(control$npsim*ncol(res$par)), control$npsim)
    }else{
      if(control$grid_type == "LHS"){
        simPoints <- lhsDesign(control$npsim, ncol(res$par))$design
      }else{
        if(control$grid_type == "grid2d"){
          simPoints <- as.matrix(expand.grid(seq(0, 1,length.out = control$npsim),
                                             seq(0, 1,length.out = control$npsim)))
        }
      }
    }
    
    if(!is.null(lower) &  !is.null(upper)){
      # rescale
      simPoints <- scale(simPoints, center = FALSE, scale = 1/(upper - lower))
      simPoints <- scale(simPoints, center = -lower, scale = FALSE)
    }
    
    if(!is.null(res$history)){
      Simu_f1 <- simulate(res$history$model[[1]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                          checkNames = FALSE, nugget.sim = 10^-8)
      Simu_f2 <- simulate(res$history$model[[2]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                          checkNames = FALSE, nugget.sim = 10^-8)
      CPF <- CPF(Simu_f1, Simu_f2, res$history$y, f1lim = control$f1lim, f2lim = control$f2lim)
    }else{
      Simu_f1 <- simulate(res$lastmodel[[1]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                          checkNames = FALSE, nugget.sim = 10^-8)
      Simu_f2 <- simulate(res$lastmodel[[2]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                          checkNames = FALSE, nugget.sim = 10^-8)
      CPF <- CPF(Simu_f1, Simu_f2, cbind(res$lastmodel[[1]]@y, res$lastmodel[[2]]@y),
                 f1lim = control$f1lim, f2lim = control$f2lim)
    }
    if(control$printVD) cat("Vorob'ev deviation: ", CPF$VD, "\n")
    
    if(add){
      points(CPF$PF, col = control$PF.points.col, pch = control$PF.pch)
      plotParetoEmp(CPF$PF, col = control$PF.line.col)
      plotParetoEmp(CPF$VE, col = control$VE.line.col, lty = 2, lwd = 2)
    }else{
      plotSymDevFun(CPF)
    }
    
    
  }
  
  if(UQ_PS){
    if(!is.null(res$history)){
      plot_uncertainty(res$history$model, lower = lower, upper = upper, resolution = control$resolution,
                       nintegpoints = control$nintegpoints, option = control$option)
    }else{
      plot_uncertainty(res$lastmodel, lower = lower, upper = upper, resolution = control$resolution,
                       nintegpoints = control$nintegpoints, option = control$option)
    }
  }
  
  if(UQ_dens){
    # reuse of conditional simulations to compute CPS if Simu_f1 is already created
    if(!exists("Simu_f1")){
      if(!is.null(res$history)){
        Simu_f1 <- simulate(res$history$model[[1]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                            checkNames = FALSE, nugget.sim = 10^-8)
        Simu_f2 <- simulate(res$history$model[[2]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                            checkNames = FALSE, nugget.sim = 10^-8)
      }else{
        Simu_f1 <- simulate(res$lastmodel[[1]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                            checkNames = FALSE, nugget.sim = 10^-8)
        Simu_f2 <- simulate(res$lastmodel[[2]], nsim = control$nsim, newdata = simPoints, cond = TRUE,
                            checkNames = FALSE, nugget.sim = 10^-8)
      }
    }
    
    non_dom_set <- NULL
    for (i in 1:control$nsim){
      non_dom_order <- nds_rank(rbind(Simu_f1[i,], Simu_f2[i,]))
      non_dom_set <- rbind(simPoints[which(non_dom_order == 1),], non_dom_set)
    }
    
    if(!is.null(res$history)){
      estDens <- ParetoSetDensity(model = res$history$model, lower = lower, upper = upper,
                                  CPS = non_dom_set)
    }else{
      estDens <- ParetoSetDensity(model = res$lastmodel, lower = lower, upper = upper,
                                  CPS = non_dom_set)
    }
    plot(estDens, display = control$displaytype)
  }
}



parplotPF_nd <- function(PF, color = NULL, add = T, lty = NULL, xlab = "", ylab = ""){
  n <- ncol(PF)
  
  if(!add){
    plot(seq(1,n), rep(0,n), pch = 3, ylim = c(min(PF)-0.05, max(PF)), xlab = xlab,
         ylab = ylab, xaxt = "n")
    axis(1, at = seq(1,n), las = 0)
    for(i in 1:n){
      abline(v = i)
    }
  }
  
  
  for(i in 1:nrow(PF)){
    if(is.null(color)){
      col <- i
    }else{
      col <- color[i]
    }
    
    if(is.null(lty)){
      ltype <- floor(i/8)+1 # 8 is the size of the standard color palette
    }else{
      ltype <- lty[i]
    }
    lines(seq(1:n), PF[i,], col = col, lty = ltype)
  }
}


