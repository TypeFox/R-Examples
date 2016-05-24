################################################################################
##
## $Id: matchedPortfolio.R 404 2007-04-19 16:27:55Z enos $
##
## Basic matchedPortfolio class.
##
## This matchedPortfolio class represents an "original" portfolioBasic that
## has been matched with one or more portfolios according to a specified
## formula.
##
## For simplicity I've included methods on matchedPortfolioCollection
## here as well, but should really be separated.
##
################################################################################

setMethod("show",
          signature(object = "matchedPortfolio"),
          function(object){
            
            ## displays essential information

            plural <- ifelse(ncol(object@matches) > 1, "s", "")

            cat("An object of class \"matchedPortfolio\":\n\n",
                "slots:\n\n", sep = "")
            
            cat("formula: "); print(object@formula); cat("\n\n")

            cat("original:\n",
                sprintf("  %i %s", nrow(object@original@weights),
                        "positions."), "\n",
                sprintf("  %.3f %s", performance(object@original)@ret,
                        object@original@ret.var), "\n\n"

                )
            
            cat(paste("matches:", sep = ""), "\n",
                sprintf("  %i %s",
                        nrow(object@matches),
                        "positions."), "\n",
                sprintf("  %i %s%s", ncol(object@matches),
                        "matched portfolio", plural), "\n",
                sprintf("  %.3f %s", totalReturn(object),
                        object@original@ret.var), "\n\n"
                )

            cat("Difference in returns:", performance(object@original)@ret
                - totalReturn(object), "\n")
            
          }
)

setMethod("summary",
          signature(object = "matchedPortfolio"),
          function(object){

            ## Really need to refactor this .orig.perf function
            ## business.
            
            orig.perf <- performance(object@original)@ret
            all.matched.perf <- totalReturn(object, output = "individual")
            matched.perf <- mean(all.matched.perf, na.rm = TRUE)
            perf.ranks <- rank(c(orig.perf, all.matched.perf[!is.na(all.matched.perf)]))
            percentile <- perf.ranks[1] / length(perf.ranks)
            
            cat("Matched portfolio summary: ", object@original@name, "\n\n",
                ncol(object@matches), " matches using ", object@method, " matching.\n\n",
                "Matched portfolio returns:\n",
                sep = "")

            print(summary(all.matched.perf))

            cat("\n",
                "Original portfolio return: ", format(orig.perf),
                ", with ", performance(object@original)@missing.return, " NAs.\n",
                "Original return relative to matches: ", format(orig.perf - matched.perf), "\n\n",
                "Original portfolio outperformed ", format(round(100 * percentile)), "% of matches.\n",
                sep = ""
                )
          }
)


## Do I really need this method?

setMethod("totalReturn",
          signature(object = "matchedPortfolio"),
          function(object, output = "pooled"){

            ## "output" is the format for returning performance.  The
            ## options are "pooled" and "individual".  "pooled"
            ## returns a scalar, the mean of the returns for all
            ## matched portfolios.  "individual" returns a vector of
            ## the returns for each portfolio

            stopifnot(validObject(object))
            stopifnot(output %in% c("individual", "pooled"))

            res <-
              sapply(data.frame(object@matches),
                     function(x){
                       sum(object@original@weights$weight[match(row.names(object@matches), object@original@weights$id)] *
                           object@original@data[[object@original@ret.var]][match(x, object@original@data$id)]) })
            
            
#            ## "matches" is a weight matrix.  "returns" is a vector of
#            ## returns. "output" is either "individual" or "pooled" and
#            ## specifies whether to return a vector of the returns for each
#            ## matched portfolio or the mean return

#            matches <- object@matches
#            returns <-
#              object@original@data[[object@original@ret.var]][match(rownames(matches),
#                                                                    object@original@data$id)]
#            ## calculates returns for each stock

#            res <- apply(matches, 2, function(x){x * returns})
            
#            ## calculates sum for each matched portfolio
            
#            res <- apply(res, 2, sum, na.rm = TRUE)
            
            if(identical(output, "pooled")){
              res <- mean(res, na.rm = TRUE)
            }

            invisible(res)

          }
          )

setMethod("performance",
          signature(object = "matchedPortfolio"),
          function(object){
            invisible(totalReturn(object, output = "individual"))
          }
          )

setMethod("contribution",
          signature(object = "matchedPortfolio", contrib.var = "character"),
          function(object, contrib.var, buckets = 5){
            
          }
          )

setMethod("exposure",
          signature(object = "matchedPortfolio", exp.var = "character"),
          function(object, exp.var, output = "pooled"){

            ## "exp.var" is a vector of variables on which to
            ## calculate exposure.  The "data" slot of the "original"
            ## portfolio must contain a column for every value in
            ## "exp.var".  Returns a list where each element in the
            ## list is the name of an element in "exp.var".

            stopifnot(validObject(object))

            ## Right now, only allow a single exposure variable that's
            ## numeric.

            data <- object@original@data
            
            ## Cannot calculate performance if there are no "exp.var"
            ## values.  "data" must contain information for at least 1
            ## stock.

            stopifnot(
                      length(exp.var) > 0,
                      nrow(data) > 0,
                      all(exp.var %in% names(data))
                      )

            exp.list <- list()

            ## Construct weight vector for matches:

            weight.matches <- object@original@weights$weight[match(row.names(object@matches),
                                                                   object@original@weights$id)]
            
            ## calculates exposures for all elements of "exp.var"

            for(ev in exp.var){

              ## 2007-04-17: This is screwed up, because the handling
              ## of numeric and character/factor exposures should be
              ## the same for different output types.  As I look at
              ## this now, there is no notion of "pooled" for numeric
              ## exposures; furthermore, use of colMeans here was
              ## incorrect and produced wrong result.  Now I use
              ## colSums correctly, but multiple results will always
              ## be returned for matched portfolios with more than one
              ## match.
              
              if(is.numeric(data[[ev]])){ ## numeric variables
                
                exp <- apply(object@matches, 2, function(x){
                  weight.matches * data[[ev]][match(x, data$id)]})

                exp.list[[ev]] <- colSums(exp)
                
              }
              else if(is.character(data[[ev]]) || is.factor(data[[ev]])){ # character variables

                ## builds a matrix to store the results of each call
                ## to "tapply".  "matrix" contains as many rows as
                ## their are matched portfolios and as many columns as
                ## there are levels of the current exposure variable,
                ## "ev"
                
                exp <- matrix(0,
                              nrow = dim(object@matches)[[2]],
                              ncol = length(levels(data[[ev]])),
                              dimnames = list(1:dim(object@matches)[[2]],
                                levels(data[[ev]])[order(levels(data[[ev]]))])
                              )

                for(i in 1:ncol(object@matches)){

                  exp[i,] <- tapply(weight.matches, data[[ev]][match(object@matches[,i], data$id)], sum,
                                    simplify = TRUE)

                  ## Make NA's 0 to make this old code work.
                  
                  exp[i,is.na(exp[i,])] <- 0

                }
                
                if(identical(output, "pooled")){
                  exp <- colMeans(exp)
                }
                
                exp.list[[ev]] <- exp

              }
            }

            exp.list

          }
          )

setMethod("plot",
          signature(x = "matchedPortfolio", y = "missing"),
          function(x, type = "returns"){

            data <- .matched.perf.df(x)

            
            panel <- function(x, subscripts, ...){
              
              panel.histogram(x, ...)
              data.sub <- data[subscripts, , drop = FALSE]
              stopifnot(sum(!data.sub$match) == 1)
              orig.ret <- data.sub$return[!data.sub$match]
              percentile <- round(100 * rank(data.sub$return)[!data.sub$match] / max(rank(data.sub$return)))
              
              grid.segments(x0 = unit(orig.ret, "native"), y0 = unit(0, "native"),
                            x1 = unit(orig.ret, "native"), y1 = unit(1, "npc"),
                            gp = gpar(lty = "dashed"))
              grid.text(paste(format(percentile), "% ", sep = ""),
                        x = unit(orig.ret, "native"), y = unit(1, "npc") - unit(1, "lines"),
                        just = "right", gp = gpar(cex = 0.8))
#              grid.text(paste(" ", format(100 - percentile), "%", sep = ""),
#                        x = unit(orig.ret, "native"), y = unit(1, "npc") - unit(1, "lines"),
#                        just = "left", gp = gpar(cex = 0.8))

            }

            main <- paste(x@original@name, "vs. matched portfolios")
            histogram( ~ return, data = data, panel = panel, nint = 50, type = "percent", main = main)
          }
          )

################################################################################
##
## Helper/utility functions
##
################################################################################

## .matched.perf.df
##
## Helper function that collects performance and summary info for a
## matchedPortfolio object.  elements of the collection

.matched.perf.df <- function(object){

  stopifnot(is(object, "matchedPortfolio"))
  
  exp.orig <- sum(object@original@weights$weight, na.rm = TRUE)
  num.orig <- sum(!is.na(object@original@weights$weight))

  ## All the matches should have the same total net exposure, and
  ## should only differ from the net exposure of the original
  ## portfolio by omitted positions. Only break if there don't appear
  ## to be any matches at all.
  
  stopifnot(nrow(object@matches) > 0)
  exp.match <- sum(object@original@weights$weight[match(row.names(object@matches),
                                                        object@original@weights$id)])
  num.match <- nrow(object@matches)
  
  data <- rbind(data.frame(match = FALSE, return = performance(object@original)@ret,
                           omitted.treatment = object@omitted.treatment,
                           omitted.control = object@omitted.control,
                           exposure = exp.orig, num = num.orig
                           ),
                
                data.frame(match = TRUE,  return = totalReturn(object, output = "individual"),
                           omitted.treatment = NA, omitted.control = NA,
                           exposure = exp.match, num = num.match

                           ))
  invisible(data)
}

.match.as.portfolioBasic <- function(object, match.num){
  
  stopifnot(is(object, "matchedPortfolio"),
            match.num <= ncol(object@matches))

  weights <- object@original@weights
  weights$id <- object@matches[,match.num][match(weights$id,
                                                 row.names(object@matches))]
  row.names(weights) <- weights$id
  
  p <- new("portfolioBasic", data = object@original@data, symbol.var = "symbol")
  p@weights <- weights

  p
  
}
