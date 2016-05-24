## Yang Lu Yang.Lu@williams.edu

## plot methods for brinson class

setMethod("plot",
          signature(x = "brinson", y = "missing"),
          function(x,
                   var = "sector",
                   type = "exposure"){
            
            stopifnot(type %in% c("exposure", "return"))            
            ## type == "exposure"
            if (type == "exposure"){
              temp <- exposure(x, var = var)
              len.row <- nrow(temp)
              df <- data.frame(Name = rep(rownames(temp), 2),
                               Value = c(temp[ , 1], temp[ , 2]),
                               Type = c(rep("Portfolio", len.row),
                                 rep("Benchmark", len.row)))
              if (class(x@universe[[var]])[1] != "numeric"){
                df <- within(df,
                             Name <- factor(Name, levels = rev(levels(x@universe[[var]]))))} else {
                               df <- within(df,
                                            Name <- factor(Name, levels = unique(as.factor(df$Name))))}
              
              ## plot
              temp.plot <- .bar.plot(df = df,
                                     type = "Exposure",
                                     title = "Exposure -- Portfolio vs. Benchmark")
              print(temp.plot)
              
            } else {
              ## type == "returns"
              df <- data.frame(Name = rep(names(x@ret.port), 2),
                               Value = c(x@ret.port, x@ret.bench),
                               Type = c(rep("Portfolio", length(x@ret.port)),
                                 rep("Benchmark", length(x@ret.bench))))

              if (class(x@universe[[var]])[1] != "numeric"){
                df <- within(df,
                             Name <- factor(Name, levels = rev(levels(x@universe[[var]]))))}else {
                               df <- within(df,
                                            Name <- factor(Name, levels = unique(as.factor(df$Name))))}
      
              
              temp.plot <- .bar.plot(df = df,
                                     type = "Return",
                                     title = "Return -- Portfolio vs. Benchmark")
              print(temp.plot)
            }
          }
          
          )

## plot methods for brinsonMulti class

setMethod("plot",
          signature(x = "brinsonMulti", y = "missing"),
          function(x,
                   var = "sector",
                   type = "exposure"){
            
            stopifnot(type %in% c("exposure", "return"))
            if (type == "exposure"){
              .df <- list()
              dates <- x@date.var
              len <- length(dates)
              temp <- exposure(x, var = var)
              len.var <- nrow(temp[[1]])
              
              for (k in 1:len){
                .df[[k]] <- data.frame(Date = rep(as.character(dates[k]), 2 * len.var),
                                       Name = rep(rownames(temp[[1]], 2)),
                                       Value = c(temp[[1]][ , k], temp[[2]][ , k]),
                                       Type = c(rep("Portfolio", len.var),
                                         rep("Benchmark", len.var)))
              }
              .df <- do.call(rbind,.df)

              if (class(x@universe[[1]]@universe[[var]])[1] != "numeric"){
                .df <- within(.df,
                              Name <- factor(Name, levels = rev(levels(x@universe[[1]]@universe[[var]]))))}else {
                                .df <- within(.df,
                                              Name <- factor(Name, levels = unique(as.factor(.df$Name))))}
              
              ## plot
              temp.plot <- .facet.plot(df = .df,
                                       type = "Exposure",
                                       title = "Exposure across Periods")
              print(temp.plot)
              
            } else {
              ## type == "returns"
              .df <- list()
              dates <- x@date.var
              len <- length(dates)
              len.var <- nrow(x@ret.port)
              for (k in 1:len){
                .df[[k]] <- data.frame(Date = rep(as.character(dates[k]), len.var),
                                       Name = rep(rownames(x@ret.port), 2),
                                       Value = c(x@ret.port[, k], x@ret.bench[, k]),
                                       Type = c(rep("Portfolio", len.var),
                                         rep("Benchmark", len.var)))
              }
              
              .df <- do.call(rbind,.df)

              if (class(x@universe[[1]]@universe[[var]])[1] != "numeric"){
                .df <- within(.df,
                              Name <- factor(Name, levels = rev(levels(x@universe[[1]]@universe[[var]]))))}else {
                                .df <- within(.df,
                                              Name <- factor(Name, levels = unique(as.factor(.df$Name))))}

              ## plot

              temp.plot <- .facet.plot(df = .df,
                                       type = "Return",
                                       title = "Return across Periods")
              print(temp.plot)
            }
          }
          )

## plot for regression class object

setMethod("plot",
          signature(x = "regression", y = "missing"),
          function(x,
                   var = "sector",
                   type = "exposure"){

            stopifnot(type %in% c("exposure", "return"))
            
            if (type == "exposure"){
              temp <- exposure(x, var = var)
              len.row <- nrow(temp)
              df <- data.frame(Name = rep(rownames(temp), 2),
                               Value = c(temp[ , 1], temp[ , 2]),
                               Type = c(rep("Portfolio", len.row),
                                 rep("Benchmark", len.row)))
              
              if (class(x@universe[[var]])[1] != "numeric"){
                df <- within(df,
                             Name <- factor(Name, levels = rev(levels(x@universe[[var]]))))} else {
                               df <- within(df,
                                            Name <- factor(Name, levels = unique(as.factor(df$Name))))}
              

              ## plot
              temp.plot <- .bar.plot(df = df,
                                     type = "Exposure",
                                     title = "Exposure -- Portfolio vs. Benchmark")
              print(temp.plot)
              
            } else {
              
              ## type == "return", no graph for a single-period
              ## regression attribution.
              print("No graph for single-period regression attribution.")
              returns(x)
            }
          }
          )


## plot for regressionMulti class object
setMethod("plot",
          signature(x = "regressionMulti", y = "missing"),
          function(x,
                   var = "sector",
                   type = "exposure"){

            stopifnot(type %in% c("exposure", "return"))
            if (type == "exposure"){
              .df <- list()
              dates <- x@date.var
              len <- length(dates)
              temp <- exposure(x, var = var)
              len.var <- nrow(temp[[1]])

              Date <- NULL
              Value <- NULL
              Type <- NULL
              rm(Date, Value, Type)
              
              for (k in 1:len){
                .df[[k]] <- data.frame(Date = rep(as.character(dates[k]), 2 * len.var),
                                       Name = rep(rownames(temp[[1]], 2)),
                                       Value = c(temp[[1]][ , k], temp[[2]][ , k]),
                                       Type = c(rep("Portfolio", len.var),
                                         rep("Benchmark", len.var)))
              }
              .df <- do.call(rbind,.df)

              
              if (class(x@universe[[1]]@universe[[var]])[1] != "numeric"){
                .df <- within(.df,
                              Name <- factor(Name, levels = rev(levels(x@universe[[1]]@universe[[var]]))))}else {
                                .df <- within(.df,
                                              Name <- factor(Name, levels = unique(as.factor(.df$Name))))}
              
              ## plot
              temp.plot <- .facet.plot(df = .df,
                                       type = "Exposure",
                                       title = "Exposure across periods")
              print(temp.plot)
              
            } else {
              ## type == "returns"
              len.row <- length(x@date.var)
              df <- data.frame(Date = rep(as.character(x@date.var), 2),
                               Value = c(x@portfolio.ret[1, ],
                                 x@benchmark.ret[1, ]),
                               Type = c(rep("Portfolio", len.row),
                                 rep("Benchmark", len.row)))
              df[1:12, 2] <- cumprod(df[1:12, 2] + 1) - 1
              df[13:24, 2] <- cumprod(df[13:24, 2] + 1) - 1

              temp.plot <- ggplot(df, aes(x = Date, y = Value, col = Type, group = Type)) +
                geom_line(aes(linetype = Type)) +
                  scale_y_continuous()+ geom_hline(yintercept = 0) + 
                    ylab("Return") + xlab("Date") +
                      ggtitle("Portfolio Performance") + 
                      theme(panel.background = element_blank(),
                            ## title = "Portfolio Performance",
                            axis.line = element_blank(), ## theme_blank(),
                            panel.grid.minor = element_blank(),
                            panel.grid.major = element_blank(),
                            axis.text.x = element_text(angle = 90, hjust = 0.5),
                            plot.background = element_rect(fill = NA, colour = NA))
                              
              print(temp.plot)
              
            }
          }
          
          )

