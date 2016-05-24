## Yang Lu Yang.Lu@williams.edu

## added in brinson attribution based on the input category
## (e.g. sectors) on Jan 1, 2013

## returns method for brinson class

setMethod("returns",
          signature(object = "brinson"),
          function(object,
                   ...){

            ## round to certain digits
            options(digits = 3)
            
            ## returns by category
            portret <- object@ret.port
            benchret <- object@ret.bench
            portwt <- object@weight.port
            benchwt <- object@weight.bench
            
            cat.allocation <- (portwt - benchwt) * benchret
            cat.selection <- (portret - benchret) * benchwt
            cat.interaction <- (portret - benchret) * (portwt - benchwt)
            
            cat.ret <- cbind(cat.allocation,
                             cat.selection,
                             cat.interaction)
            colnames(cat.ret) <- c("Allocation", "Selection", "Interaction")
            cat.ret <- rbind(cat.ret, apply(cat.ret, 2, sum))
            rownames(cat.ret)[nrow(cat.ret)] <- "Total"
            
            ## overall brinson attribution
            q1 <- object@q1
            q2 <- object@q2
            q3 <- object@q3
            q4 <- object@q4

            asset.allocation <- q2 - q1
            stock.selection <- q3 - q1
            interaction <- q4 - q3 - q2 + q1
            active.ret <- q4 - q1
            
            ret.mat <- matrix(NA, nrow = 4, ncol = 1)
            ret.mat[1, 1] <- asset.allocation
            ret.mat[2, 1] <- stock.selection
            ret.mat[3, 1] <- interaction
            ret.mat[4, 1] <- active.ret
            
            colnames(ret.mat) <- as.character(unique(object@universe[[object@date.var]]))
            rownames(ret.mat) <- c("Allocation Effect",
                                   "Selection Effect",
                                   "Interaction Effect",
                                   "Active Return")

            ## organize output
            output.list <- list()
            output.list[[1]] <- cat.ret * 10000
            output.list[[2]] <- ret.mat
            names(output.list) <- c("Attribution by category in bps",
                                    "Aggregate")
            
            return(output.list)
          }
          )

## returns method for brinsonMulti class with three different types of
## compounding - arithmetic, geometric and linking coefficient
## approach.

setMethod("returns",
          signature(object = "brinsonMulti"),
          function(object,
                   type = "geometric",
                   ...){

            stopifnot(type %in% c("arithmetic", "linking", "geometric"))
                      
            ## three types - default is geometric, the other two are
            ## arithmetic and linking
            
            active.return <- object@brinson.mat[1,] - object@brinson.mat[4,]
            allocation <- object@brinson.mat[3,] - object@brinson.mat[4,]
            selection <- object@brinson.mat[2,] - object@brinson.mat[4,]
            interaction <- active.return - allocation - selection
            
            ari.raw <- rbind(allocation, selection, interaction, active.return)
            rownames(ari.raw) <- c("Allocation", "Selection",
                                   "Interaction", "Active Return")
            
            if (type == "arithmetic"){
              ari.agg <- .aggregate(object, ari.raw)
              ari.list <- .combine(ari.raw, ari.agg)
              return(ari.list)
            }

            if (type == "linking"){
              port.ret.overtime <- apply(object@brinson.mat + 1, 1,prod)[1]
              bench.ret.overtime <- apply(object@brinson.mat + 1, 1, prod)[4]
              act.return <- port.ret.overtime - bench.ret.overtime

              T <- dim(object@brinson.mat)[2]

              A.natural.scaling <- act.return / T /
                ((port.ret.overtime) ^ (1 / T) - (bench.ret.overtime) ^ (1 / T))

              names(A.natural.scaling) <- NULL

              C <- (act.return -
                    A.natural.scaling * sum(object@brinson.mat[1,] -
                                            object@brinson.mat[4,])) /
                (sum((object@brinson.mat[1,] - object@brinson.mat[4,])^2))
              alpha <- C * (object@brinson.mat[1,] - object@brinson.mat[4,])

              B.linking <- A.natural.scaling + alpha
              
              linking.raw <- t(sapply(1:4, function(i){ari.raw[i,] * B.linking})) 
              rownames(linking.raw) <- c("Allocation", "Selection",
                                         "Interaction", "Active Return")
              linking.agg <- .aggregate(object, linking.raw)
              linking.list <- .combine(linking.raw, linking.agg)
              return(linking.list)
            }

            if (type == "geometric"){
              temp.mat <- apply(object@brinson.mat + 1, 1, prod) - 1
              allocation <- temp.mat[3] - temp.mat[4] ## q2 - q1
              selection <- temp.mat[2] - temp.mat[4] ## q3 - q1
              
              ## q4 - q3 - q2 + q1
              interaction <- temp.mat[1] - temp.mat[2] - temp.mat[3] + temp.mat[4]
              
              active.ret <- temp.mat[1] - temp.mat[4] ## q4 - q1
              
              ret.mat <- matrix(NA, nrow = 4, ncol = 1)
              ret.mat[1, 1] <- allocation
              ret.mat[2, 1] <- selection
              ret.mat[3, 1] <- interaction
              ret.mat[4, 1] <- active.ret
              
              colnames(ret.mat) <- paste(c(min(unique(as.character(object@date.var))),
                                           max(unique(as.character(object@date.var)))),
                                         collapse = ", ")
              rownames(ret.mat) <- c("Allocation", "Selection",
                                     "Interaction", "Active Return")
              geo.list <- .combine(ari.raw, ret.mat)
              return(geo.list)
            }
            
          }
          )

## returns for regression class object
setMethod("returns",
          signature(object = "regression"),
          function(object,
                   ...){
            
            ## round to certain digits
            options(digits = 3)
            no.row <- length(object@reg.var)
            ret.mat <- matrix(NA, nrow = no.row + 4, ncol = 1)

            j <- 1
            for (i in 1:no.row){
              col.name <- object@universe[[object@reg.var[i]]]
              if (class(col.name)[1] != "numeric"){
                ret.mat[i, 1] <- sum(object@contrib[j:(j + length(levels(col.name)) - 1)])
                j <- j + length(levels(col.name))
              } else {
                ret.mat[i, 1] <- object@contrib[j]
                j <- j + 1
              }
            }
                      
            ret.mat[no.row + 1, 1] <- object@act.ret - sum(object@contrib)
            ret.mat[no.row + 2, 1] <- object@portfolio.ret
            ret.mat[no.row + 3, 1] <- object@benchmark.ret
            ret.mat[no.row + 4, 1] <- object@act.ret

            colnames(ret.mat) <- as.character(unique(object@universe[[object@date.var]]))
            rownames(ret.mat) <- c(object@reg.var,
                                   "Residual",
                                   "Portfolio Return",
                                   "Benchmark Return",
                                   "Active Return")
            return(ret.mat)
          }
          )


## returns for regressionMulti class object
setMethod("returns",
          signature(object = "regressionMulti"),
          function(object,
                   type = "geometric",
                   ...){
            
            ## 3 types - default is geometric, the other two are
            ## arithmetic and linking.
            stopifnot(type %in% c("arithmetic", "linking", "geometric"))
            
            ## raw attribution
            raw <- sapply(1:length(object@date.var), function(i){returns(object@universe[[i]])})
            rownames(raw) <- rownames(returns(object@universe[[1]]))
            colnames(raw) <- unique(object@date.var)

            if (type == "arithmetic"){
              agg <- matrix(apply(raw, 1, sum))
              colnames(agg) <- paste(c(min(unique(as.character(object@date.var))),
                                       max(unique(as.character(object@date.var)))),
                                     collapse = ", ")
              rownames(agg) <- rownames(returns(object@universe[[1]]))
              ari.list <- .combine(raw, agg)
              return(ari.list)
            }
            
            if (type == "linking"){
              T <- length(object@date.var)
              portfolio.ret <- apply(object@portfolio.ret + 1, 1, prod)
              benchmark.ret <- apply(object@benchmark.ret + 1, 1, prod)
              act.ret <- portfolio.ret - benchmark.ret
              A <- act.ret / T / (portfolio.ret ^ (1 / T) -
                                  benchmark.ret ^ (1 / T))
              C <- (act.ret - sum(object@act.ret * A)) / sum(object@act.ret^2)
              alpha <- C * object@act.ret
              B.linking <- A + alpha
              .mat <- returns(object)[["Raw"]]
              .no.row <- nrow(.mat)
              linking.raw <- t(sapply(1:.no.row, function(i){.mat[i,] * B.linking})) 
              rownames(linking.raw) <- rownames(.mat)
              colnames(linking.raw) <- colnames(.mat)
              linking.raw <- linking.raw[c(-.no.row + 1, -.no.row + 2), ]
              
              linking.agg <- matrix(apply(linking.raw, 1, sum))
              colnames(linking.agg) <- paste(c(min(unique(as.character(object@date.var))),
                                       max(unique(as.character(object@date.var)))),
                                     collapse = ", ")
              rownames(linking.agg) <- rownames(linking.raw)
              linking.list <- .combine(linking.raw, linking.agg)
              return(linking.list)
            }
            
            if (type == "geometric"){
              agg <- matrix(apply(raw + 1, 1, prod) - 1)
              .no.row <- nrow(agg)
              agg[.no.row, 1] <- agg[.no.row - 2, 1] - agg[.no.row - 1, 1]
              rownames(agg) <- rownames(returns(object@universe[[1]]))
              colnames(agg) <- paste(c(min(unique(as.character(object@date.var))),
                                       max(unique(as.character(object@date.var)))),
                                     collapse = ", ")
              geo.list <- .combine(raw, agg)
              return(geo.list)
            }
          }
          )

