#' Mixed Membership Visualization
#' 
#' 
#' \code{vizTheta} plots \eqn{\theta}, the parameters which govern the sub-population distributions of variables in a mixed membership model.
#'  The parameter \eqn{\theta_{j,k}} specifies the distribution of variable j for complete members of sub-population k. The estimated
#' parameters from the given model are shown in black in each plot; the parameters from a comparison model (if available) are shown in red.
#' Each row of plots represents a single variable, and each column of the plots represents a sub-population. 
#' 
#' This is the function called by the plot generic function for \code{mixedMemModel} objects.  
#'  
#' @param model the \code{mixedMemModel} object that will be plotted.
#' @param compare an array of the same dimensions as model$theta which contains values to compare against the fitted value. 
#' @param main the main figure title.
#' @param varNames a vector of strings specifying labels for each variable.
#' @param groupNames a vector of strings specifying labels for each sub-population.
#' @param nrow the number of rows in each plot. If the argument is not specified, all variables will appear in one plot.
#' @param indices a vector which indicates specific variables to plot. If the argument is not specified, all variables will be plotted. If the number of variables
#' to plot is greater than nrow, then multiple plots will be produced.
#' @param fitNames the names of the models plotted.
vizTheta = function(model, compare = NULL, main = "Estimated Theta",
                    varNames = NULL, groupNames = NULL,nrow = NULL, fitNames = NULL, indices = NULL) {
  
  # Internal Variables to set
  h.space <- .25
  v.space <- .1
  
  if(is.null(varNames)) {
    varNames <- paste("Var", c(1:model$J))
  }
  
  if(is.null(groupNames)) {
    groupNames <- paste("Group", c(1:model$K))
  }
  
  if (is.null(nrow)) {
    nrow <- model$J
  }
  
  if(is.null(indices)) {
    indices <- c(1:model$J)
  }
  
  if(is.null(fitNames)) {
    if(is.null(compare)) { 
      fitNames <- paste("Model", 1)
    } else {
      fitNames <- paste("Model", c(1:2))
    }
  }
  
  graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, model$K), mar = rep(.1,4))
  count = 1
  for(j in indices)
  {
    if(model$dist[j]=="multinomial"|model$dist[j]=="rank")
    {
      for(k in 1:model$K)
      {
        graphics::plot(model$theta[j,k,], type = "p", lwd = 2, col = "black", ylim = c(-v.space, 1+v.space), xlim = c(h.space, model$Vj[j] + h.space),
             yaxt = "n", xaxt = "n", pch = 16)
        if(!is.null(compare)) {
          graphics::points(c(1:model$Vj[j]), compare[j,k,c(1:model$Vj[j])], col = "red", pch = 4, lwd = 1.5)
        }
        if (k == 1) {
          graphics::mtext(varNames[j], line = 3, side = 2, cex = .7)
          graphics::axis(side = 2, at = c(0,.5,1), labels = c(0,.5,1))
        }
        if(count == indices[length(indices)]| (count %% nrow) == 0) {
          graphics::mtext(paste(groupNames[k], sep = " "), line = .2, side = 1, cex = 1-min(model$J,10)*.4)
        }
      }
    } else if (model$dist[j] == "bernoulli") {
      for(k in 1:model$K) {
        graphics::plot(model$theta[j,k,], type = "p", lwd = 2, col = "black", ylim = c(-v.space,1+ v.space), xlim = c(h.space, 1 + h.space),
             yaxt = "n", xaxt = "n", pch = 16)
        if(!is.null(compare)) {
          graphics::points(c(1:model$Vj[j]), compare[j,k,], col = "red", pch = 4, lwd = 1.5)
        }
        if(k == 1){
          graphics::mtext(varNames[j], line = 3, side = 2, cex = .7)
          graphics::axis(side = 2, at = c(0,.5,1), labels = c(0,.5,1))
        }
        if(count == indices[length(indices)]| (count %% nrow) == 0) {
          graphics::mtext(paste(groupNames[k], sep = " "), line = .2, side = 1, cex = .8)
        }
      }
    } 
    if((count %% nrow) == 0) {
      graphics::title(main = main, outer = T, cex = 1.2)      
      graphics::par(fig = c(0, 1, 0, 1), oma = c(0,5,0,1), mar = rep(0, 4), new = T)
      
      graphics::plot(0, 0, type = "n", bty = "n", xaxt ="n", yaxt = "n")
      
      if(is.null(compare)){
        graphics::legend("bottom", legend = fitNames, pch = 19,
               col = "black", cex = .8)
        
      } else {
        graphics::legend("bottom", legend = fitNames,
               pch = c(19, 4), col = c("black", "red"), ncol = 2, cex = .8)
      }
      graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, model$K), mar = rep(.1,4))
    }
    count = count + 1
  }
  graphics::title(main = main, outer = T, cex = 1.2)      
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0,5,0,1), mar = rep(0, 4), new = T)
  
  graphics::plot(0, 0, type = "n", bty = "n", xaxt ="n", yaxt = "n")
  
  if(is.null(compare)){
    graphics::legend("bottom", legend = fitNames, pch = 19,
           col = "black", cex = .8)
    
  } else {
    graphics::legend("bottom", legend = fitNames,
           pch = c(19, 4), col = c("black", "red"), ncol = 2, cex = .8)
  }
  graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, model$K), mar = rep(.1,4))
}


#' Mixed Membership Visualization
#' 
#' 
#' \code{vizMem} plots estimates for the group membership scores of each individual. This is the function called by the \code{mixedMemModel} class
#' generic plot function. 
#' 
#' The estimates plotted are the normalized \eqn{\phi}, which are the posterior means from the variational distribution. The
#' estimated group membership scores are  shown in black, and the estimates from a comparison model (if available) are shown in red. Each plot
#' represents an individual. 
#' 
#' This is the function called by the plot generic function for mixedMemModel objects.
#' 
#' @param model the \code{mixedMemModel} object that will be plotted.
#' @param compare a matrix of the same dimensions as \code{model$phi} which contains parameters to compare
#'  against the fitted model. 
#' @param main the main figure title.
#' @param groupNames a vector specifying labels for each sub-population.
#' @param nrow the number of rows in each plot.
#' @param ncol the number of columns in each plot.
#' @param indices the specific individuals which will be shown in the plot. If the argument is left blank, all individuals will be plotted.
#' If the number of individuals to be plotted is larger than \code{nrow * ncol} then multiple plots will be produced.
#' @param fitNames a vector of length 2 containing strings which correspond to the names of the models (fitted and comparison).
#' @export
vizMem <- function(model, compare = NULL, main = "Estimated Membership",
                   nrow = NULL, ncol = NULL,
                   indices = NULL, groupNames = NULL, fitNames = NULL) {
  # Internal Variables
  pch.list = c(0:14)
  h.space = .25
  v.space = .1
  
  if (is.null(indices)) {
    indices <- c(1:model$Total)   
  }
  
  if (is.null(nrow)) {
    nrow = floor(sqrt(length(indices)))
  }
  
  if(is.null(ncol)) {
    ncol = ceiling(length(indices)/nrow)
  }
  
  if (is.null(groupNames)) {
    groupNames = paste("Group", c(1:model$K))
  }
  
  if(is.null(fitNames)) {
    if(is.null(compare)) { 
      fitNames <- paste("Model", 1)
    } else {
      fitNames <- paste("Model", c(1:2))
    }
  }
  
  mem.est <- model$phi / rowSums(model$phi)
  
  graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, ncol), mar = rep(.1,4))
  count <- 0
  for(i in indices)
  {
    count <- count + 1
    
    graphics::plot(mem.est[i,], type = "p", lwd = 2, col = "black", ylim = c(-v.space,1 + v.space), xlim = c(h.space, model$K + h.space),
         yaxt = "n", xaxt = "n", pch = pch.list)
    graphics::text(0, 1, labels = i , cex = .9, adj = c(0, 1), pos = 4)
    if((count %% ncol) == 1) {
      graphics::axis(2, at = c(0,.5,1), labels = c(0,.5,1), cex = 1)
    }
    
    if(!is.null(compare)){
      graphics::points(c(1:model$K), compare[i,], col = "red", pch = pch.list)
    }
    if((count %%(nrow*ncol)) == 0)
    {
      graphics::title(main = main, outer = T, cex = 1.2)
      graphics::mtext("Group Membership", side = 2, outer = T, line = 3)
      graphics::mtext("Groups", side = 1, outer = T, line = 0)
      
      graphics::par(fig = c(0, 1, 0, 1), oma = c(0,5,0,1), mar = rep(0, 4), new = T)
      
      graphics::plot(0, 0, type = "n", bty = "n", xaxt ="n", yaxt = "n")
      if(is.null(compare)){
        graphics::legend("bottom", legend = paste(fitNames[1], groupNames), pch = pch.list[1:model$K],
               col = rep("black", model$K), cex = .8, ncol = model$K)
        
      } else {
        graphics::legend("bottom", legend = c(paste(fitNames[1], groupNames),paste(fitNames[2], groupNames) ),
               pch = pch.list[1:model$K], col = rep(c("black", "red"), each = model$K), ncol = model$K, cex = .8)
      }
      graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, ncol), mar = rep(.1,4))
    }
  }
  graphics::title(main = main, outer = T, cex = 1.2)
  graphics::mtext("Group Membership", side = 2, outer = T, line = 3)
  graphics::mtext("Groups", side = 1, outer = T, line = 0)
  
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0,5,0,1), mar = rep(0, 4), new = T)
  
  graphics::plot(0, 0, type = "n", bty = "n", xaxt ="n", yaxt = "n")
  
  if(is.null(compare)){
    graphics::legend("bottom", legend = paste(fitNames[1], groupNames), pch = pch.list[1:model$K],
           col = rep("black", model$K), cex = .8, ncol = model$K)
    
  } else {
    graphics::legend("bottom", legend = c(paste(fitNames[1], groupNames),paste(fitNames[2], groupNames) ),
           pch = pch.list[1:model$K], col = rep(c("black", "red"), each = model$K), ncol = model$K, cex = .8)
  }
  graphics::par(oma = c(3,5,3,1), mfrow = c(nrow, ncol), mar = rep(.1,4))
}
