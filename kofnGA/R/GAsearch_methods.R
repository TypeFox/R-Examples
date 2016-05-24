summary.GAsearch <- function(object, ...) {
  
  out <- list()
  out$ng <- length(object$old$obj) - 1  #-Number of generations.
  out$nu <- dim(unique(object$pop))[1]  #-Number of unique solutions in the final population.
  out$bs <- object$bestsol  #-Best solution.
  out$bsg <- which(object$old$obj == min(object$old$obj))[1] - 1  #-Gen at which best was found.
  out$OFvals <- cbind(average = object$old$avg[c(1, out$ng + 1)], 
                      minimum = object$old$obj[c(1, out$ng + 1)])
  rownames(out$OFvals) <- c("Initial population", "Final population")
  
  class(out) <- "summary.GAsearch"
  out
  
}


print.summary.GAsearch <- function(x, ...) {
  
  cat("Genetic algorithm search,", x$ng, "generations\n")
  cat("Number of unique solutions in the final population:", x$nu, "\n\n")
  cat("Objective function values:\n")
  print(x$OFvals)
  cat("\nBest solution (found at generation ", x$bsg, "):\n", sep="")
  cat(x$bs, "\n")
  
}


print.GAsearch <- function(x, ...) {
  
  ng <- length(x$old$obj) - 1  #-Number of generations.
  ps <- dim(x$pop)[1]  #-Population size.
  k <- dim(x$pop)[2]  #-Subset size.
  
  cat("Genetic algorithm search,", ng, "generations\n\n")
  cat("$pop:", "       ", ps, "x", k, " matrix\n", sep="")
  cat("$obj:", "       ", ps, "x1  vector\n", sep="")
  cat("$old$obj:", "   ", ng + 1, "x1  vector\n", sep="")
  cat("$old$avg:", "   ", ng + 1, "x1  vector\n", sep="")
  cat("$old$best:", "  ", ng + 1, "x1  vector\n\n", sep="")
  
  cat("$bestsol:\n", x$bestsol, "\n\n")
  cat("$bestobj:\n", x$bestobj, "\n\n")
  
}


plot.GAsearch <- function(x, type=c("l","l"), lty=c(1,1), pch=c(-1,-1), col=c("blue", "red"), lwd=c(1,1), ...) { 
  
  ng <- length(x$old$obj) - 1
  lb = min(c(x$old$obj, x$old$avg))
  ub = max(c(x$old$obj, x$old$avg))
  bestloc = which(x$old$obj == min(x$old$obj))[1] - 1
  
  matplot(0:ng, cbind(x$old$avg, x$old$obj), type=type, lty=lty, pch=pch, col=col, lwd=lwd, 
          xlab="Generation", ylab="Objective function value", ...)
  
  txt <- c("Pop average", "Pop best", "Last improvement")   
  # Need to handle plot inputs as either numeric or character.
  if(is.character(lty)) leglty <- c(lty, "blank") else leglty <- c(lty, 0)
  if(is.character(pch)) legpch <- c(pch, "x") else legpch <- c(pch, 19)
  if(is.character(col)) legcol <- c(col, "black") else legcol <- c(col,1)
  
  points(bestloc, x$old$obj[bestloc + 1], pch=legpch[3])
  
  legend("topright", legend=txt, pch=legpch, col=legcol,
         lwd=c(lwd,1), bty="n", lty=leglty, cex=0.8)  
  
}