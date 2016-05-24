## show method
## Yang Lu Yang.Lu@williams.edu

## The show method displays the first part of the summary, which
## contains only the essential info about the brinson model.

setMethod("show",
          signature(object = "brinson"),
          function(object){
            
            cat(sprintf("Period:                              %-s",
                        paste(unique(as.character(object@universe[[object@date.var]])),
                              collapse = ", ")), "\n",
                sprintf("Methodology:                         %-s",
                        paste("Brinson")), "\n",
                sprintf("Securities in the portfolio:         %-s",
                        length(which(object@universe[[object@portfolio.weight]] > 0))), "\n",
                sprintf("Securities in the benchmark:         %-s",
                        length(which(object@universe[[object@bench.weight]] > 0))), "\n",
                sep = ""
                )
            
            cat("\n")
          }
          )


## The show method displays the first part of the summary, which
## contains only the essential info about the brinson model.

setMethod("show",
          signature(object = "brinsonMulti"),
          function(object){
            
            cat(sprintf("Period starts:                       %-s",
                        paste(min(unique(as.character(object@date.var))),
                              collapse = ", ")), "\n",
                sprintf("Period ends:                         %-s",
                        paste(max(unique(as.character(object@date.var))),
                              collapse = ", ")), "\n",
                sprintf("Methodology:                         %-s",
                        paste("Brinson")), "\n",
                sprintf("Securities in the portfolio:         %-s",
                        do.call(mean, lapply(1:length(object@date.var),
                                             function(i){length(which(object@universe[[i]]@universe[[object@portfolio.weight]] > 0))}))), "\n",
                sprintf("Securities in the benchmark:         %-s",
                        do.call(mean, lapply(1:length(object@date.var),
                                             function(i){length(which(object@universe[[i]]@universe[[object@bench.weight]] > 0))}))), "\n",
                sep = ""
                )
            
            cat("\n")
          }
          )


## show method for regression class object

setMethod("show",
          signature(object = "regression"),
          function(object){
            
            cat(sprintf("Period:                              %-s",
                        paste(unique(as.character(object@universe[[object@date.var]])),
                              collapse = ", ")), "\n",
                sprintf("Methodology:                         %-s",
                        paste("Regression")), "\n", 
                sprintf("Securities in the portfolio:         %-s",
                        length(which(object@universe[[object@portfolio.weight]] > 0))), "\n",
                sprintf("Securities in the benchmark:         %-s",
                        length(which(object@universe[[object@benchmark.weight]] > 0))), "\n",
                sep = ""
                )
            
            cat("\n")

          }
          )


## show method for regressionMulti class object

setMethod("show",
          signature(object = "regressionMulti"),
          function(object){
            
            cat(sprintf("Period starts:                       %-s",
                        paste(min(unique(as.character(object@date.var))),
                              collapse = ", ")), "\n",
                sprintf("Period ends:                         %-s",
                        paste(max(unique(as.character(object@date.var))),
                              collapse = ", ")), "\n",
                sprintf("Methodology:                         %-s",
                        paste("Regression")), "\n",
                sprintf("Securities in the portfolio:         %-s",
                        do.call(mean, lapply(1:length(object@date.var),
                                             function(i){length(which(object@universe[[i]]@universe[[object@portfolio.weight]] > 0))}))), "\n",
                sprintf("Securities in the benchmark:         %-s",
                        do.call(mean, lapply(1:length(object@date.var),
                                             function(i){length(which(object@universe[[i]]@universe[[object@benchmark.weight]] > 0))}))), "\n",
                sep = ""
                )
            
            cat("\n")

          }
          )

