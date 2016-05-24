
#setClass("bbo",representation(minCost='list', bestHabitat="list"))

#setClass("Paths", representation(best.paths='list', path.scores="list"))



summary.bbo <- function(object, ...){
  digits <- max(5, getOption('digits') - 2)
 
  cat(" ::summary of BBO run:: ",
      "-----------",
      "Properties:",
      "-----------", 
      sep = "\n")

  cat("numVar: ", object$prop$numVar,
      "\npopSize: ", object$prop$popSize,
      "\nmaxGen: ", object$prop$maxGen,
      "\nKeep: ", object$prop$KEEP,
      "\npMutate: ", object$prop$pMutate,
      "\npModify: ", object$prop$pModify,
      "\norderDep: ", object$prop$orderDep, "\n")

  cat("\n--------------",
      "Best Solution:",
      "--------------",
      "Best habitat:", 
      sep = "\n")

  print(round(object$minCost$bestMember, digits))
  cat("Best value/minimum cost:", 
       sep = "\n")
  print(round(object$minCost$bestValue, digits))
  cat("\n------------",
     "Generations:",
     "------------",
     "Average population value for each generation:",
     sep = "\n")
  print(round(object$bestHabitat$itersAvg, digits))
  cat("Best habitat for each generation:", 
      sep = "\n")
  print(round(object$bestHabitat$itersBestMember, digits))
  cat("Best(minimum) function cost for each generation:", 
      sep = "\n")
  print(round(object$bestHabitat$itersBestValue, digits))

  invisible(object)
}

plot.bbo <- function (x, plot.type = c("itersAvg", "itersBestValue"), ...) {

 if (identical(plot.type[1], "itersBestValue")) {
    #plot.new()
    if( length(grep("R version 3", R.Version()$version.string)) ){
 	plot(x$bestHabitat$itersBestValue, type = "b", xlab = "#generation", ylab = "Best cost among the population", main = "BBO: Best cost among population members at each generation", col = 'blue', cex.main = 1.5, cex = 0.8)
    }else{
	plot(x$bestHabitat$itersBestValue, type = "b", xlab = "#generation", ylab = "Best cost among the population", main = "BBO: Best cost among population members at each generation", col = 'blue', cex.main = 0.9, cex = 0.8)
    }
 
 }else if (identical(plot.type[1], "itersAvg")) {
    #plot.new()
    if( length(grep("R version 3", R.Version()$version.string)) ){
    	plot(x$bestHabitat$itersAvg, type = "b", xlab = "#generation", ylab = "Average cost of population members", main = "BBO: Population average at each generation", col = "blue", cex.main = 1.5, cex = 0.8)
    }else{
	plot(x$bestHabitat$itersAvg, type = "b", xlab = "#generation", ylab = "Average cost of population members", main = "BBO: Population average at each generation", col = "blue", cex.main = 0.9, cex = 0.8)
    }


 }else {
   warning("'plot.type' does not correspond to any plotting type", immediate. = TRUE)
 }

}
