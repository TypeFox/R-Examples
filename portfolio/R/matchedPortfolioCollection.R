################################################################################
##
## $Id: matchedPortfolioCollection.R 389 2007-01-10 04:28:44Z enos $
##
## Collection (list) of 'matchedPortfolio' objects.
##
################################################################################

setMethod("summary",
          signature(object = "matchedPortfolioCollection"),
          function(object){

            data <- .matched.coll.perf.df(object)

            data.matched <- data[data$match,]
            data.orig    <- data[!data$match,]

            data.summary <- merge(data.orig[c("period","return","omitted.treatment","omitted.control","exposure","num")],
                                  aggregate(data.matched[c("return","exposure","num")], by = list(period = data.matched$period), mean),
                                  by = "period")
            
            names(data.summary) <- c("period","orig","omit.trt","omit.ctrl","orig.exp","orig.num","match","match.exp","match.num")
            data.summary$orig.sum <- cumsum(data.summary$orig)
            data.summary$match.sum <- cumsum(data.summary$match)
            data.summary$diff <- data.summary$orig - data.summary$match
            
            data.pct <- lapply(split(data, data$period),
                               function(x){rank(x$return)[!x$match] / max(rank(x$return)) })
            data.percentile <- data.frame(period = names(data.pct), percentile = round(100 * unlist(data.pct)))

            data.summary <- merge(data.summary, data.percentile, by = "period")

            cat("Matched portfolio collection summary ", "\n\n", sep = "")
            print(data.summary)

          }
)

setMethod("plot",
          signature(x = "matchedPortfolioCollection", y = "missing"),
          function(x, type = "returns"){

            ## Wouldn't it be nice if the plot function or any
            ## (single) function returned a grob that we could use?
            ## Much of this code is duplicated for 'matchedPortfolio'...
            
            data <- .matched.coll.perf.df(x)

            ## data <- data[order(data$period, decreasing = TRUE),]
            ## data$period <- as.factor(data$period)
            ## data$period <- factor(data$period, levels = rev(levels(data$period)))

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

            ## Just take first name for simplicity.
            
            main <- paste(x@data[[1]]@original@name, "vs. matched portfolios")
            histogram( ~ return | period, data = data, panel = panel, nint = 100, main = main)
          }
          )


setMethod("matching",
          signature(object = "data.frame"),
          function(object,
                   universe   = NULL,
                   covariates = NULL,
                   method     = "greedy",
                   n.matches  = 1,
                   exact      = character(),
                   by.var,
                   in.var,
                   weight.var = NULL,
                   ret.var,
                   ...
                   ){

            stopifnot(!is.null(by.var) &&
                      !is.null(in.var))
            
            matchedPortfolio.from.df <- function(df, ...){

              p.tomatch <- new("portfolioBasic",
                               data       = df,
                               type       = "relative",
                               in.var     = in.var,
                               weight.var = weight.var,
                               ret.var    = ret.var,
                               size       = "all", ...)
              m <- matching(p.tomatch,
                            universe   = universe,
                            covariates = covariates,
                            method     = method,
                            n.matches  = n.matches,
                            exact      = exact)
            }
            
            res <- lapply(split(object, object[[by.var]]), function(x) { matchedPortfolio.from.df(x, ...) })
            m <- new("matchedPortfolioCollection", data = res)
            invisible(m)

          }
          )

## .matched.coll.perf.df
##
## Helper function that collects performance and summary info for all
## elements of the collection

.matched.coll.perf.df <- function(object){

  stopifnot(is(object, "matchedPortfolioCollection"))
  
  data <- NULL
  
  for(n in names(object@data)){
    data.sub <- .matched.perf.df(object@data[[n]])
    data.sub$period <- n
    data <- rbind(data.sub, data)
  }

  data <- data[order(data$period),]
  
  invisible(data)
}

