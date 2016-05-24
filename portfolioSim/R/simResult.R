################################################################################
##
## $Id: simResult.R 401 2007-04-19 16:18:51Z enos $
##
## Methods defined for class 'simResult'.
##
################################################################################

setMethod("summary",
          signature(object = "simResult"),
          function(object, start.period = NULL, end.period = NULL){
            
            cat("\nSimulation summary:\n\n")
            
            ## The simResult object contains a list of single period
            ## results.  For now, collect the pertinent information
            ## for this summary into a data frame.

            if("basic" %in% object@type){
              
              ret <- perfSummaryDf(object, start.period = start.period, end.period = end.period)
              
              ## Haven't dealt with NA returns yet.
              
              stopifnot(all(!is.na(ret$ret)))
              
              ## The 'ret' data frame really ought be sorted at this
              ## point, but since we're relying on sorting in
              ## calculations like drawdowns, sort again just in case.

              ret <- ret[order(ret$period),]
            
              if(nrow(ret) > 0){

                periods <- ret$period
                periods <- periods[order(periods)]

                ## Don't call min and max since we don't require that
                ## these functions are defined on the class of the
                ## period data.
              
                period.min <- periods[1]
                period.max <- periods[length(periods)]

                ## Some convenience functions for dealing with choosing
                ## between percent and bps displays.
              
                pct.thresh <- 0.01

                .ret.scale <- function(x){
                  stopifnot(length(x) == 1)
                  ifelse(abs(x) > pct.thresh, x * 100, x * 100 * 100)
                }
                .ret.tag <- function(x){
                  stopifnot(length(x) == 1)
                  ifelse(abs(x) > pct.thresh, "%", "bps")
                }
                
                .commify <- function(x){
                  stopifnot(length(x) == 1)
                  format(round(x), big.mark = ",")
                }
                
                ## Summary statistics.
                
                ## total.ret     <- prod(ret$ret + 1) - 1
                total.ret     <- sum(ret$ret)
                total.profit  <- sum(ret$profit)
                
                mean.ret      <- mean(ret$ret)
                ann.mean.ret  <- mean.ret * object@freq
                
                vol.ret       <- sd(ret$ret)
                ann.vol.ret   <- vol.ret * sqrt(object@freq)
                
                sharpe        <- ann.mean.ret / ann.vol.ret
                
                ## Find best and worst periods.
                
                best.period.df   <- ret[which.max(ret$ret),]
                worst.period.df  <- ret[which.min(ret$ret),]
                
                ## Find worst drawdown.
                
                ret$nav  <- cumprod(1 + ret$ret)
                dd <- 1 - ret$nav / cummax(ret$nav)
                dd.max <- max(dd)
                
                ## Do a little extra work so that we can use periods
                ## that don't have comparison operators.
                
                ## dd.start really represents the end-of-period moment
                ## in time when the drawdown occurs.  It does _not_
                ## include the return for that period.
              
                dd.end.idx   <- max(which.max(dd), length(dd) - which.max(rev(dd)) +1)
                dd.end       <- ret$period[dd.end.idx]
                dd.start.idx <- max(which(dd == 0 & seq(dd) < dd.end.idx))
                dd.start     <- ret$period[dd.start.idx]
                dd.ret       <- sum(ret$ret[dd.start.idx:dd.end.idx])
                dd.profit    <- sum(ret$profit[dd.start.idx:dd.end.idx])
                dd.range.tag <- ifelse(dd.start.idx == dd.end.idx,
                                       format(dd.end),
                                       sprintf("%s to %s", format(dd.start), format(dd.end)))

                ## Turnover
                
                mean.size          <- mean(ret$size)
                mean.equity        <- mean(ret$equity)
                mean.gross.equity  <- mean(ret$gross.equity)
                universe.turnover  <- sum(ret$universe.turnover)
                turnover           <- sum(ret$turnover)
                mean.turnover      <- mean(ret$turnover)
                ann.mean.turnover  <- mean.turnover * object@freq
              
                holding.period.months <- ((mean.gross.equity * 2) / ann.mean.turnover) * 12
                
                ## Missingness
                
                mean.missing.price <- mean(ret$missing.price)
                mean.missing.return <- mean(ret$missing.return)
                
                
                cat("\n")
                
                cat(sprintf("Start:   %s", format(period.min)), "\n",
                    sprintf("End:     %s", format(period.max)), "\n\n",
                    
                    sprintf("                                    profit   return"), "\n",
                    sprintf("---------------------------------------------------"), "\n",
                    sprintf("Total:                   %17s   %6.1f %-3s",
                            .commify(total.profit),
                            .ret.scale(total.ret),
                            .ret.tag(total.ret)), "\n",
                    sprintf("Sharpe:                  %17s   %6.1f",
                            "",
                            sharpe), "\n",
                    sprintf("---------------------------------------------------"), "\n",
                    sprintf("Mean return:             %17s   %6.1f %-3s",
                            "",
                            .ret.scale(mean.ret),
                            .ret.tag(mean.ret)
                            ), "\n",
                    sprintf("Mean return (ann):       %17s   %6.1f %-3s",
                            "",
                            .ret.scale(ann.mean.ret),
                            .ret.tag(ann.mean.ret)
                            ), "\n",
                    sprintf("Volatility:              %17s   %6.1f %-3s",
                            "",
                            .ret.scale(vol.ret),
                            .ret.tag(vol.ret)
                            ), "\n",
                    sprintf("Volatility (ann):        %17s   %6.1f %-3s",
                            "",
                            .ret.scale(ann.vol.ret),
                            .ret.tag(ann.vol.ret)
                            ), "\n",
                    sprintf("---------------------------------------------------"), "\n",
                    sprintf("Best period:             %17s   %6.1f %-3s (%s)",
                            .commify(best.period.df$profit),
                            .ret.scale(best.period.df$ret),
                            .ret.tag(best.period.df$ret),
                            format(best.period.df$period)), "\n",
                    sprintf("Worst period:            %17s   %6.1f %-3s (%s)",
                            .commify(worst.period.df$profit),
                            .ret.scale(worst.period.df$ret),
                            .ret.tag(worst.period.df$ret),
                            format(worst.period.df$period)), "\n",
                    sprintf("Worst drawdown:          %17s   %6.1f %-3s (%s)",
                            .commify(dd.profit),
                            .ret.scale(dd.ret),
                            .ret.tag(dd.ret),
                            dd.range.tag
                            ), "\n",
                    sprintf("---------------------------------------------------"), "\n",
                    sprintf("Mean size:               %17s",
                            .commify(mean.size)), "\n",
                    sprintf("Mean equity:             %17s",
                            .commify(mean.equity)), "\n",
                    sprintf("Mean gross equity:       %17s",
                            .commify(mean.gross.equity)), "\n",
                    sprintf("Universe turnover:       %17s",
                            .commify(universe.turnover)), "\n",
                    sprintf("Turnover:                %17s",
                            .commify(turnover)), "\n",
                    sprintf("Mean turnover:           %17s",
                            .commify(mean.turnover)), "\n",
                    sprintf("Mean turnover (ann):     %17s",
                            .commify(ann.mean.turnover)), "\n",
                    sprintf("Holding period (mth):    %17.1f",
                            holding.period.months), "\n",
                    sprintf("---------------------------------------------------"), "\n",
                    sprintf("Mean NA weights:         %17s",
                            .commify(mean.missing.price)
                            ), "\n",
                    sprintf("Mean NA returns:         %17s",
                            .commify(mean.missing.return)
                            ), "\n\n",
                    sep = ""
                    
                    )
              }
            }

            if("detail" %in% object@type){
              per.sec <- securityPerfSummaryDf(object, start.period = start.period, end.period = end.period)

              ## First show by profit, then contrib.

              per.sec <- per.sec[order(per.sec$profit),]
              cat("Top/bottom performers by profit:\n\n")
              print(rbind(head(per.sec[!is.na(per.sec$profit),c("id","profit"),]),
                          tail(per.sec[!is.na(per.sec$profit),c("id","profit"),])))

              per.sec <- per.sec[order(per.sec$contrib),]
              
              cat("\nTop/bottom performers by contribution (%):\n\n")
              print(rbind(head(per.sec[!is.na(per.sec$contrib),c("id","contrib"),]),
                          tail(per.sec[!is.na(per.sec$contrib),c("id","contrib"),])))
  
              cat("\n\n")
            }

            if("contributions" %in% object@type){
              contrib <- contributionSummaryDf(object, start.period = start.period, end.period = end.period)
              
              cat("Mean contributions:\n\n")
              print(contrib)
            }

            if("exposures" %in% object@type){

              ## Summary of exposures data
              
            }

            if("trades" %in% object@type){

              ## Summary of trades data
              
            }

            if(!is.null(object@summary.interface)){
              summary(object@summary.interface)
            }
            
          }
          )


setMethod("plot",
          signature(x = "simResult", y = "missing"),
          function(x, y, ...){

            dots <- unlist(c(...))

            ret <- perfSummaryDf(x)

            main <- "Portfolio simulation result"
            print(xyplot(cumprod(ret + 1) - 1 ~ period,
                         data = ret,
                         type = "l",
                         main = main))
          }
          )


## Helper function for taking a bunch of return data in period.data
## objects and turning into a data frame.  

setMethod("perfSummaryDf",
          signature = signature(object = "simResult"),
          function(object, start.period = NULL, end.period = NULL){
            ret <- do.call(rbind,
                           lapply(object@data,
                                  function(x){
                                    if((is.null(start.period) || x@period.data@period > start.period) &
                                       (is.null(end.period)   || x@period.data@period <= end.period)){
                                      data.frame(period            = x@period.data@period,
                                                 start             = x@start.data@instant,
                                                 end               = x@end.data@instant,
                                                 ret               = x@period.data@performance@ret,
                                                 profit            = x@period.data@performance@profit,
                                                 turnover          = x@period.data@turnover,
                                                 universe.turnover = x@period.data@universe.turnover,
                                                 missing.price     = x@period.data@performance@missing.price,
                                                 missing.return    = x@period.data@performance@missing.return,
                                                 equity.long       = x@start.data@equity.long,
                                                 equity.short      = x@start.data@equity.short,
                                                 size.long         = x@start.data@size.long,
                                                 size.short        = x@start.data@size.short
                                                 )
                                    }
                                    else{
                                      NULL
                                    }
                                  }
                                  )
                           )
            
            if(nrow(ret) == 0) return(NULL)
            
            ## Should I really have period.label here, instead of period?

            ## It's important to note, then, that the start instant of the first
            ## period supplied by a client application refers to the moment of
            ## initial investment.  At that instant, the NAV is 1.

            ## Tack on an initial position row to the front of the returns data
            ## frame.  The start period of the first record is used as the label
            ## -- better hope this is different than the first period!

            ret$equity <- mapply(function(long, short) { mean(c(long, abs(short))) },
                                 ret$equity.long,
                                 ret$equity.short)

            ## Also tack on a notion of gross equity.

            ret$gross.equity <- ret$equity.long + abs(ret$equity.short)

            ## Size is the total number of positions (not an average of
            ## long/short.

            ret$size <- ret$size.long + ret$size.short
            
            ret <- ret[c(1, 1:nrow(ret)),]
            ret[1,] <- data.frame(period = ret$start[2],
                                  start             = NA,
                                  end               = NA,
                                  ret               = 0,
                                  profit            = 0,
                                  turnover          = 0,
                                  universe.turnover = 0,
                                  missing.price     = 0,
                                  missing.return    = 0,
                                  equity.long       = ret$equity.long[2],
                                  equity.short      = ret$equity.short[2],
                                  size.long         = ret$size.long[2],
                                  size.short        = ret$size.short[2],
                                  equity            = ret$equity[2],
                                  gross.equity      = ret$gross.equity[2],
                                  size              = ret$size[2]
                                  )

            ## We're trying to make a sensible first period tag, but
            ## if the classes of period and start are not the same,
            ## trouble will be afoot.  Do our best to deal.

            if(class(ret$period) != class(ret$start)){
              if(is.numeric(ret$period)){
                ret$period[1] <- ret$period[2] - 1
              }
            }
            
            if(ret$period[1] == ret$period[2]){
              warning(paste("First period tag and start instant are the same.",
                            "Plots and summaries may not function properly"))
              
            }

            invisible(ret)
            
          }
          )

## Helper function for taking stock-specific performance data from the
## result objects and condensing for use in a summary method.

setMethod("securityPerfSummaryDf",
          signature = signature(object = "simResult"),
          function(object, start.period = NULL, end.period = NULL){

            ## At what point does this become intractable?  I may need
            ## to do some smart bookeeping along the way here...
            
            x <- do.call(rbind,
                         lapply(object@data,
                                function(x){
                                  if((is.null(start.period) || x@period.data@period > start.period) &
                                     (is.null(end.period)   || x@period.data@period <= end.period)  &
                                     nrow(x@period.data@performance@ret.detail) > 0){
                                    x@period.data@performance@ret.detail[c("id","contrib","profit")]
                                  }
                                  else{
                                    NULL
                                  }
                                }))
            x.agg <- aggregate(x[c("contrib","profit")], by = list(id = x$id), sum, na.rm = TRUE)

            ## Turn contributions into percents

            x.agg$contrib <- x.agg$contrib * 100

            invisible(x.agg)
          }
          )

## Helper function for taking contribution data from the result
## objects and condensing for use in a summary method.

## I should probably call this something different, since it returns a
## list of data frames this time.

setMethod("contributionSummaryDf",
          signature = signature(object = "simResult"),
          function(object){

            ## At what point does this become intractable?  I may need
            ## to do some smart bookeeping along the way here...

            contrib.vars <- unique(unlist(lapply(object@data,
                                                 function(x){
                                                   names(x@period.data@contribution@data)
                                                 })))

            all.contrib <- list()

            for(contrib.var in contrib.vars){
              x <- do.call(rbind,
                           lapply(object@data,
                                  function(x){
                                    y <- x@period.data@contribution@data[[contrib.var]][c("rank","contrib","roic")]
                                    y$date <- x@period.data@period
                                    y
                                  }))

              ## Turn into percents

              x$contrib <- x$contrib * 100
              x$roic    <- x$roic * 100
              
              all.contrib[[contrib.var]] <- x
            }

            
            all.contrib <- lapply(all.contrib,
                                  function(x){
                                    aggregate(x[c("contrib","roic")], by = list(rank = x$rank), mean)
                                  }
                                  )
            
            invisible(all.contrib)
          }
          )


setMethod("loadIn",
          signature(object = "simResult", in.loc = "character", fmt = "missing"),
          function(object, in.loc){

            ## Again, we only have a single format: directories with a
            ## period.RData file.  If you find one of those, the data
            ## in its directory can be loaded up into a single period
            ## result object.

            load(paste(in.loc, "master.RData", sep = "/"))
            object@freq <- master.freq
            object@errors <- master.errors
            object@type <- master.type
            object@summary.interface <- master.summary.interface
            
            files <- grep("period\\.RData",
                          list.files(in.loc, full.names = TRUE, recursive = TRUE), value = TRUE)
            
            in.locs <- gsub("\\/[^/]*$", "", files)

            object@data <- vector("list", length(in.locs))

            if(length(in.locs) > 0){
              for(i in 1:length(in.locs)){
                print(in.locs[i])
                sp <- new("simResultSinglePeriod")
                sp <- loadIn(sp, in.loc = in.locs[i])
                object@data[[i]] <- sp
              }
            }
            invisible(object)
          }
          )


## This method saves out all the slots of a simResult other than data
## to the top level of the out.loc directory

setMethod("saveOut",
          signature(object  = "simResult",
                    type    = "missing",
                    fmt     = "missing",
                    out.loc = "character",
                    name    = "missing",
                    verbose = "logical"),
          function(object, out.loc, verbose){

            dir.create(out.loc, showWarnings = FALSE, recursive = TRUE)
            
            master.freq <- object@freq
            master.errors <- object@errors
            master.type <- object@type
            master.summary.interface <- object@summary.interface

            save(master.freq,
                 master.errors,
                 master.type,
                 master.summary.interface,
                 file = paste(out.loc, "master.RData", sep = "/")
                 )
          }
          )
