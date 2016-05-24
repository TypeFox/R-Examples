summary.rbga <- function(object, echo=FALSE, ...) {
    # should do object type checking here
    rbga.object = object
    
    output = paste(
        "GA Settings", "\n",
        "  Type                  = ", rbga.object$type, "\n",
        "  Population size       = ", rbga.object$popSize, "\n",
        "  Number of Generations = ", rbga.object$iters, "\n",
        "  Elitism               = ", rbga.object$elitism, "\n",
        "  Mutation Chance       = ", rbga.object$mutationChance, "\n",
        "\n",
        "Search Domain", "\n", sep="")
    for (var in 1:length(rbga.object$stringMin)) {
        minVar = rbga.object$stringMin[var]
        maxVar = rbga.object$stringMax[var]
        output = paste(output,
            "  Var ", var, " = [", minVar, ",", maxVar, "]\n", sep="");
    }
    output = paste(output, "\n", sep="");
    if (!is.null(rbga.object$suggestions)) {
        optionPart = paste("Suggestions", "\n", sep="");
        for (suggestion in 1:dim(rbga.object$suggestions)[1]) {
            optionPart = paste(optionPart,
                "  ", suggestion, " = (", sep="");
            for (var in 1:(dim(rbga.object$suggestions)[2]-1)) {
                optionPart = paste(optionPart,
                    rbga.object$suggestions[suggestion,var], ", ",
                    sep="");
            }
            # and the last one
            optionPart = paste(optionPart,
                rbga.object$suggestions[suggestion,dim(rbga.object$suggestions)[2]],
                ")\n", sep="");
        }
        output = paste(output, optionPart, "\n", sep="");
    }
    if (!is.null(rbga.object$population)) {
        optionPart = paste("GA Results", "\n", "  Best Solution : ", sep="");
        # ok, deal with the situation that more than one object is best
        filter = rbga.object$evaluations == min(rbga.object$evaluations);
        bestObjectCount = sum(rep(1, rbga.object$popSize)[filter]);
        if (bestObjectCount > 1) {
            bestSolution = rbga.object$population[filter,][1,];
        } else {
            bestSolution = rbga.object$population[filter,];
        }
        for (var in 1:length(bestSolution)) {
            optionPart = paste(optionPart,
                bestSolution[var], " ",
                sep="");
        }
        output = paste(output, optionPart, "\n", sep="");
    }
    if (echo) cat(output);
    invisible(output);
}

