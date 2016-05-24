setMethod(f = "auc"
          , signature(object = "bincap")
          , function(object){
            n0 <- sum((object@response
                       == object@true) * 1 )
            n1 <- sum((object@response
                       != object@true) * 1 )
            s0 <- sum(rank(object@predicted
			   , ties.method = "average"
			   , na.last = TRUE)
                      * (object@response
                         == object@true)
                      )
            return(
                   (s0 - n0 * (n0 + 1) / 2) / ( n0 * n1 )
                   )
          }
          )
setMethod(f = "auc"
          , signature(object = "multcap")
          , function(object){
            return(mean(
                        combn(levels(object@response), 2,
                              function(levels
                                       , response
                                       , predicted){
                                df <- as.data.frame(predicted) ## factor and matrix -> need data.frame
                                df$obs <- response
                                dfs <- subset(df, get("obs") %in% levels)
                                t <- levels[1]
                                aij <- auc(new("bincap"
                                                , response = factor(dfs[,"obs"]) ## to drop non-ocurring levels
                                                , predicted = dfs[,t]
                                                , true = t)
                                            )
                                t <- levels[2]
                                aji <- auc(new("bincap"
                                                , response = factor(dfs[,"obs"]) ##  to drop non-ocurring levels
                                                , predicted = dfs[,t]
                                                , true = t)
                                            )
                                Aij <- mean(c(aij,aji))
                                return(Aij)
                              }
                              , response = object@response
                              , predicted = object@predicted
                              ), na.rm = TRUE
                        )
                   )

          }
          )

