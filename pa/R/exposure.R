## Yang Lu Yang.Lu@williams.edu

## exposure method for brinson class
## post: matrix

setMethod("exposure",
          signature(object = "brinson"),
          function(object,
                   var = "sector",
                   ...){
            
            ## round to certain digits
            options(digits = 3)
                      
            ## decide whether it's categorical or continuous (split
            ## into 5 quantiles)
            if (class(object@universe[[var]])[1] != "numeric"){
              ## categorical
              expo.mat <- cbind(tapply(object@universe[[object@portfolio.weight]],
                                       object@universe[[var]],
                                       sum),
                                tapply(object@universe[[object@bench.weight]],
                                       object@universe[[var]],
                                       sum))
              
              ## updated -- show difference betw portfolio and benchmark
              expo.mat <- cbind(expo.mat, expo.mat[,1] - expo.mat[,2])
              colnames(expo.mat) <- c("Portfolio", "Benchmark", "Diff")

              return(expo.mat)
            } else {

              ## continuous (5 quantiles)
              temp <- object@universe
              temp$q <- as.factor(ceiling(rank(temp[[var]]) / nrow(temp) * 5))
              expo.mat <- cbind(tapply(temp[[object@portfolio.weight]],
                                       temp$q,
                                       sum),
                                tapply(temp[[object@bench.weight]],
                                       temp$q,
                                       sum))

              expo.mat <- cbind(expo.mat, expo.mat[,1] - expo.mat[,2])
              
              colnames(expo.mat) <- c("Portfolio", "Benchmark", "Diff")
              rownames(expo.mat)[1] <- "Low"
              rownames(expo.mat)[5] <- "High"
              return(expo.mat)
            }
          }
          )

## exposure method for brinsonMulti class
## post: list

setMethod("exposure",
          signature(object = "brinsonMulti"),
          function(object,
                   var = "sector",
                   ...){
            
            ## round to certain digits
            options(digits = 3)

            if (class(object@universe[[1]]@universe[[var]])[1] != "numeric"){
              ## categorical
              expo.list <- list()
              no.date <- length(object@date.var)
              port.mat <- sapply(1:no.date,
                                 function(i){tapply(object@universe[[i]]@universe[[object@portfolio.weight]],
                                                    object@universe[[i]]@universe[[var]],
                                                    sum)})
              colnames(port.mat) <- object@date.var
              expo.list[[1]] <- port.mat
              
              bench.mat <- sapply(1:no.date, function(i){tapply(object@universe[[i]]@universe[[object@bench.weight]],
                                                          object@universe[[i]]@universe[[var]],
                                                          sum)})
              colnames(bench.mat) <- object@date.var
              expo.list[[2]] <- bench.mat
              expo.list[[3]] <- expo.list[[1]] - expo.list[[2]]
              
              names(expo.list) <- c("Portfolio", "Benchmark", "Diff")
              return(expo.list)
            } else {
              ## continous (5 quantiles)

              expo.list <- list()
              no.date <- length(object@date.var)

              for (i in 1:no.date){
                object@universe[[i]]@universe$q <- as.factor(ceiling(rank(object@universe[[i]]@universe[[var]]) / nrow(object@universe[[i]]@universe) * 5))
              }

              port.mat <- sapply(1:no.date,
                                 function(i){tapply(object@universe[[i]]@universe[[object@portfolio.weight]],
                                                    object@universe[[i]]@universe$q,
                                                    sum)})
              colnames(port.mat) <- object@date.var
              rownames(port.mat)[1] <- "Low"
              rownames(port.mat)[5] <- "High"
              expo.list[[1]] <- port.mat
              
              bench.mat <- sapply(1:no.date,
                                  function(i){tapply(object@universe[[i]]@universe[[object@bench.weight]],
                                                     object@universe[[i]]@universe$q,
                                                     sum)})
              
              
              colnames(bench.mat) <- object@date.var
              rownames(bench.mat)[1] <- "Low"
              rownames(bench.mat)[5] <- "High"
              expo.list[[2]] <- bench.mat
              expo.list[[3]] <- expo.list[[1]] - expo.list[[2]]
              
              names(expo.list) <- c("Portfolio", "Benchmark", "Diff")
              return(expo.list)
            }
          }
          
          )


## exposure based on regression class object

setMethod("exposure",
          signature(object = "regression"),
          function(object,
                   var = NULL,
                   ...){

            ## round to certain digits
            options(digits = 3)
            stopifnot(var %in% names(object@universe))
            
            if (!is.null(var)){
            
            ## decide whether it's categorical or continuous (split
            ## into 5 quantiles)
            if (class(object@universe[[var]])[1] != "numeric"){
              ## categorical
              expo.mat <- cbind(tapply(object@universe[[object@portfolio.weight]],
                                       object@universe[[var]],
                                       sum),
                                tapply(object@universe[[object@benchmark.weight]],
                                       object@universe[[var]],
                                       sum))
              
              ## updated -- show difference betw portfolio and benchmark
              expo.mat <- cbind(expo.mat, expo.mat[,1] - expo.mat[,2])
              colnames(expo.mat) <- c("Portfolio", "Benchmark", "Diff")

              return(expo.mat)
            } else {

              ## continuous (5 quantiles)
              temp <- object@universe
              temp$q <- as.factor(ceiling(rank(temp[[var]]) / nrow(temp) * 5))
              expo.mat <- cbind(tapply(temp[[object@portfolio.weight]],
                                       temp$q,
                                       sum),
                                tapply(temp[[object@benchmark.weight]],
                                       temp$q,
                                       sum))
              ## updated -- show difference betw portfolio and benchmark
              expo.mat <- cbind(expo.mat, expo.mat[,1] - expo.mat[,2])
              colnames(expo.mat) <- c("Portfolio", "Benchmark", "Diff")

              rownames(expo.mat)[1] <- "Low"
              rownames(expo.mat)[5] <- "High"
              return(expo.mat)
            }
          } else {
            return("Please enter a variable name")
          }
          }
          )



## regressionMulti exposure

setMethod("exposure",
          signature(object = "regressionMulti"),
          function(object,
                   var = NULL,
                   ...){
            
            ## round to certain digits
            options(digits = 3)
            stopifnot(var %in% names(object@universe[[1]]@universe))

            if (!is.null(var)){
            if (class(object@universe[[1]]@universe[[var]])[1] != "numeric"){
              ## categorical
              expo.list <- list()
              no.date <- length(object@date.var)
              port.mat <- sapply(1:no.date,
                                 function(i){tapply(object@universe[[i]]@universe[[object@portfolio.weight]],
                                                    object@universe[[i]]@universe[[var]],
                                                    sum)})
              colnames(port.mat) <- object@date.var
              expo.list[[1]] <- port.mat
              
              bench.mat <- sapply(1:no.date, function(i){tapply(object@universe[[i]]@universe[[object@benchmark.weight]],
                                                          object@universe[[i]]@universe[[var]],
                                                          sum)})
              colnames(bench.mat) <- object@date.var
              expo.list[[2]] <- bench.mat
              expo.list[[3]] <- expo.list[[1]] - expo.list[[2]]
              
              names(expo.list) <- c("Portfolio", "Benchmark", "Diff")

              return(expo.list)
            } else {
              ## continous (5 quantiles)

              expo.list <- list()
              no.date <- length(object@date.var)

              for (i in 1:no.date){
                object@universe[[i]]@universe$q <- as.factor(ceiling(rank(object@universe[[i]]@universe[[var]]) / nrow(object@universe[[i]]@universe) * 5))
              }

              port.mat <- sapply(1:no.date,
                                 function(i){tapply(object@universe[[i]]@universe[[object@portfolio.weight]],
                                                    object@universe[[i]]@universe$q,
                                                    sum)})
              colnames(port.mat) <- object@date.var
              rownames(port.mat)[1] <- "Low"
              rownames(port.mat)[5] <- "High"
              expo.list[[1]] <- port.mat
              
              bench.mat <- sapply(1:no.date,
                                  function(i){tapply(object@universe[[i]]@universe[[object@benchmark.weight]],
                                                     object@universe[[i]]@universe$q,
                                                     sum)})
              
              
              colnames(bench.mat) <- object@date.var
              rownames(bench.mat)[1] <- "Low"
              rownames(bench.mat)[5] <- "High"
              expo.list[[2]] <- bench.mat
              expo.list[[3]] <- expo.list[[1]] - expo.list[[2]]
              
              names(expo.list) <- c("Portfolio", "Benchmark", "Diff")

              return(expo.list)
            }
          } else {
            return("Please enter a variable name")
          }
          }
          )
