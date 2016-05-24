setMethod("show",
          signature(object = "movG"),
          function(object)
          {
            sumRes <- list(FunCall= object@FunCall[[1]],
                           maxValues = object@maxValues,
                           meanNumValues =  mean(object@nValues),
                           nObs = length(object@observedPhe),
                           nNAvalues = sum(is.na(object@observedPhe)),
                           correlation = object@correlation,
                           regCoef = coef(object@adjModel)[2],
                           dimFieldMap = c(nrow= max(object@row),
                             ncol = max(object@col)))
            ##
            cat("\nAdjustment by function:" ,sumRes$FunCall, "\n")
            cat("\nMaximum possible number of values in the grid:", sumRes$maxValues,"\n")
            cat("\nMean number of values in grid:", round(sumRes$meanNumValues,2),"\n")
            cat("\nNumber of observation:", sumRes$nObs,", number of NA-observations:",
                sumRes$nNA,"\n")
            cat("\nCoefficient of correlation between moving means\nand observed phe. values:",
                round(sumRes$correlation,2),"\n")
            cat("\nRegression coefficient of moving means regressed\non observed phe. values :",
                round(sumRes$regCoef,2),"\n")
            cat("\nExperimental layout:\n",
                "\trows:", sumRes$dimField[1],"\n\tcolumns:",sumRes$dimField[2],"\n" )
          }
          )


