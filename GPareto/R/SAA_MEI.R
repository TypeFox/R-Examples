## ' Sample Average Approximation for several multi-objective Expected Improvement with respect to the
## ' current Pareto front. To avoid numerical instabilities, the new point is evaluated only if it is not too close to an existing observation.
## ' 
## ' @title SAA MO Expected Improvement with m objectives
## ' 
## ' @param x a vector representing the input for which one wishes to calculate EHI,
## ' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
## ' @param critcontrol List with six arguments: 
## '        \code{nb_samp} number of random samples from the posterior distribution;
## '        \code{seed} seed used for the random samples;
## '        \code{type} choice of multi-objective improvement function: "maximin" or "hypervolume";
## '        \code{refPoint} (optional) reference point for Hypervolume Expected Improvement,
## '        Options for the \code{\link[GPareto]{checkPredict}} function: \code{threshold} and \code{distance} are available to avoid numerical issues occuring when adding points too close to the existing ones.
## ' @param type "SK" or "UK" (by default), depending whether uncertainty related to trend estimation 
## '        has to be taken into account. 
## ' @param paretoFront matrix corresponding to the Pareto Front (one output per column).
## ' @return The SAA hypervolume improvement at \bold{x}.
## ' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim from which \code{EHI} is an extension.
## ' @export
## ' @importFrom MASS mvrnorm
## ' @references Svenson, J. D., & Santner, T. J. (2010). Multiobjective Optimization of Expensive Black-Box
## ' Functions via Expected Maximin Improvement. Technical Report, 
## ' 
## ' @examples
## ' #---------------------------------------------------------------------------
## ' # SAAmEI surface associated with the "P1" problem at a 15 points design
## ' #---------------------------------------------------------------------------
## ' \donttest{
## ' set.seed(25468)
## ' library(DiceDesign)
## ' 
## ' n_var <- 2 
## ' f_name <- "P1" 
## ' n.grid <- 101
## ' test.grid <- expand.grid(seq(0, 1,, n.grid), seq(0, 1,, n.grid))
## ' n_appr <- 15 
## ' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
## ' response.grid <- t(apply(design.grid, 1, f_name))
## ' Front_Pareto <- t(nondominated_points(t(response.grid)))
## ' mf1 <- km(~., design = design.grid, response = response.grid[,1])
## ' mf2 <- km(~., design = design.grid, response = response.grid[,2])
## ' 
## ' SAAmEI_grid <- apply(test.grid, 1, SAA_mEI, model = list(mf1, mf2),
## '                      critcontrol = list(nb_samp = 20, type = "maximin"))
## ' 
## ' filled.contour(seq(0, 1,, n.grid), seq(0, 1,, n.grid), matrix(SAAmEI_grid, n.grid),
## '                main = "Expected Maximin Improvement", xlab = expression(x[1]),
## '                ylab = expression(x[2]), color = terrain.colors, nlevels = 50,
## '                plot.axes = {axis(1); axis(2);
## '                             points(design.grid[,1], design.grid[,2],pch = 21,bg = "white")
## '                             }
## '               )
## ' }

SAA_mEI <- function(x, model,
                    critcontrol = list(nb.samp, seed = 42, type = "maximin", refPoint = NULL),
                    type = "UK", paretoFront = NULL){
  
  n.obj <- length(model)
  d <- model[[1]]@d
  x.new <- matrix(x, 1, d)

  if(is.null(paretoFront)){
    observations <- matrix(0, model[[1]]@n, n.obj)
    for (i in 1:n.obj) observations[,i] <- model[[i]]@y
    paretoFront <- t(nondominated_points(t(observations)))
  }
  
  refPoint <- critcontrol$refPoint
  
  nb.samp <- critcontrol$nb.samp
  if(is.null(nb.samp)){
    nb.samp <- 50
  }
  
  seed <- critcontrol$seed
  if(is.null(seed)){
    seed <- 42
  }
  
  
  Improvement <- Maximin_Improvement
  if(critcontrol$type == "hypervolume"){
    Improvement <- Hypervolume_improvement
    if (is.null(refPoint)){
      refPoint <- matrix(apply(paretoFront, 2, max) + 1, 1, n.obj) ### May be changed? 
      cat("No refPoint provided, ", signif(refPoint, 3), "used \n")
    } 
  }
  
  mu    <- rep(NaN, n.obj)
  sigma <- rep(NaN, n.obj)
  for (i in 1:n.obj){    
    pred     <- predict(object=model[[i]], newdata=x.new, type=type, checkNames = FALSE, light.return = TRUE)
    mu[i]    <- pred$mean
    sigma[i] <- pred$sd
  }
  
  ## A new x too close to the known observations could result in numerical problems
  if(checkPredict(x, model, type = type, distance = critcontrol$distance, threshold = critcontrol$threshold)){
    return(-1)
  }else{
    # Set seed to have a deterministic function to optimize
    # see http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/
    if (exists(".Random.seed", .GlobalEnv))
      oldseed <- .GlobalEnv$.Random.seed
    else
      oldseed <- NULL
    
    set.seed(seed)
    
    Samples <- mvrnorm(n = nb.samp, mu, diag(sigma))
    
    if (!is.null(oldseed)) 
      .GlobalEnv$.Random.seed <- oldseed
    else
      rm(".Random.seed", envir = .GlobalEnv)
    
    
    Res <- apply(Samples, 1, Improvement, front = paretoFront, refPoint = refPoint)
    
    return(mean(Res))
  }
  
}

Hypervolume_improvement <- function(point, front, refPoint){
  Hi <- 0
  if(!is_dominated(t(rbind(point,front)))[1]){
    Hi <- hypervolume_indicator(t(rbind(point,front)), t(front), refPoint)  
  }
  return(-Hi)
}

# ref_point added to have the same arguments than Hypervolume_improvement
Maximin_Improvement <- function(point, front, refPoint){
  Em <- 0
  if(!is_dominated(t(rbind(point,front)))[1]){
    Em <- -Inf
    for(i in 1:nrow(front)){
      tmp <- Inf
      for(j in 1:ncol(front)){
        tmp <- min(tmp, point[j] - front[i,j])
      }
      Em <- max(Em, tmp)
    }
  }
  return(-Em)
}


