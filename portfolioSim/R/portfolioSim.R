################################################################################
##
## $Id: portfolioSim.R 1322 2009-05-04 13:02:19Z enos $
##
## Methods for the main simulation object.
##
################################################################################

setMethod("initialize",
          signature(.Object = "portfolioSim"),
          function(.Object, ...){
            .Object <- callNextMethod()

            types <- .Object@out.type

            ## Type 'default'
            
            if("default" %in% types){
              types <- append(types, c("basic", "trades"))

              if(length(.Object@exp.var) > 0)
                types <- append(types, "exposures")
              if(length(.Object@contrib.var) > 0)
                types <- append(types, "contributions")

              types <- types[types != "default"]
            }

            ## Type 'all'
            
            if("all" %in% types){
              types <- append(types, c("basic", "detail", "contributions",
                                       "exposures", "trades", "portfolio"))
              types <- types[types != "all"]
            }

            ## Type 'lean'
            
            if("lean" %in% types){
              types <- append(types, c("basic"))

              if(length(.Object@exp.var) > 0)
                types <- append(types, "exposures")
              if(length(.Object@contrib.var) > 0)
                types <- append(types, "contributions")

              types <- types[types != "lean"]
            }

            types <- types[!duplicated(types)]
            
            .Object@out.type <- types

            .Object
          }
          )

setMethod("runSim",
          signature(object = "portfolioSim", verbose = "logical"),
          function(object, verbose = FALSE){

            ## Call initialize again to make sure out.type is what we
            ## want it to be.

            object <- initialize(object)

            ## Our holdings is the portfolio we are anlysing in this
            ## simulation.  Holdings at the start of each period are
            ## the same as holdings at the end of the previous period.

            ## At the start of the sim, curr.holdings gets the
            ## contents of the start.holdings slot.  Note that the
            ## data slot of curr.holdings is populated at the start of
            ## each period in the main loop.

            curr.holdings <- object@start.holdings
            
            ## Sort the periods.  The validity function for
            ## 'portfolioSim' should ensure that the following call to
            ## order works.

            periods <- object@periods[order(object@periods$period),]
            
            ## Set up the master result.  We preallocate a list large
            ## enough to store a single period result for each period.
            
            master.result <- new("simResult",
                                 freq = object@freq,
                                 data = vector("list", dim(periods)[1]),
                                 type = object@out.type
                                 )

            ## Create working copy of summary interface

            sim.summary <- object@summary.interface
            
            for(i in seq_len(dim(periods)[1])){

              prev.holdings <- curr.holdings
              
              try.result <- try({

                ## What period are we in?
                
                period <- periods$period[i]
                if(verbose) cat("Processing period:", format(period), "\n")

                ## There are two interfaces from which we gather data.
                ## This first interface supplies us with all of the
                ## data we need for the period.  Note that the result
                ## is encapsulated in its own class.

                if(verbose) cat("Getting data...")
                sim.data <- getSimData(object@data.interface, period,
                                       verbose = verbose)
                data <- sim.data@data
                if(verbose) cat("done.\n")

                ## Now that we have our period data, supply it to the
                ## holdings portfolio object.

                curr.holdings@data <- data
                curr.holdings@price.var <- "start.price"

                ## At the current level of implementation it's
                ## preferable use average price for turnover
                ## calculations.

                data$average.price <- (data$start.price + data$end.price) / 2
                
                ## We're ready at this point to clean up holdings by
                ## removing any positions that are no longer in the
                ## investable universe.  Simplification here is that I
                ## have to be able to price the securities of the
                ## positions I'm removing.  Securities I can't price
                ## stay in the portfolio and will contribute to the
                ## missing.price count.
                ##
                ## I can fail to price because I'm missing either the
                ## start or end price for the period (both are
                ## required), or the security isn't even in our data
                ## object.  The latter case I believe I'll take care
                ## of with a new missing count called missing.data.
                ## For now, both cases are caught by the check below
                ## on NA price.

                hld        <- curr.holdings@shares
                hld$price  <- data$average.price[match(hld$id, data$id)]

                hld.remove <- hld[!is.na(hld$price) &
                                  ! hld$id %in% data$id[data$universe],]

                ## The universe.turnover item contains the total
                ## currency amount of positions removed due to a
                ## changing investable universe.  There should be a
                ## penalty associated with this housecleaning.  The
                ## penalty wouldn't appear in turnover but would be
                ## accounted for in profit.

                universe.turnover <-
                  sum(abs(hld.remove$price * hld.remove$shares))
                  
                ## Finally, remove the positions from holdings and
                ## recalc weights.
                
                curr.holdings@shares <- subset(curr.holdings@shares,
                                               ! id %in% hld.remove$id)
                curr.holdings <- calcWeights(curr.holdings)
                
                ## The second interface supplies us with the set of
                ## trades we want to do for this period.  Again, the
                ## result is encapsulated in its own class.

                if(verbose) cat("Getting trades...")
                sim.trades <- getSimTrades(object@trades.interface,
                                           period,
                                           curr.holdings,
                                           sim.data,
                                           verbose = verbose)
                trades <- sim.trades@trades@trades
                if(verbose) cat("done.\n")
                
                ## We use average of start/end price in turnover
                ## calculations.
                
                trades$price <- data$average.price[match(trades$id, data$id)]

                ## At this point we "fill" the trades submitted by the
                ## trades interface for this period.  Currently, there
                ## is only one parameter to control filling behaviour,
                ## which allows fills on a fixed percentage of the
                ## day's volume.

                pct <- object@fill.volume.pct

                if(nrow(trades) > 0 && !pct %in% Inf){

                  trades$volume <- data$volume[match(trades$id, data$id)]
                  trades <- trades[!is.na(trades$volume),]
                  
                  trades$shares <- ifelse(trades$shares > trades$volume * (pct / 100),
                                          trades$volume * (pct / 100),
                                          trades$shares)

                  trades <- trades[!is.na(trades$shares) & trades$shares > 0,]
                }

                ## What should we do about trades without prices?
                ## Here, they don't make it into the sum.
                
                turnover <- sum(abs(trades$price * trades$shares), na.rm = TRUE)
                
                ## Calculate performance.  Before doing so, ensure
                ## that each stock in the portfolio at least has an
                ## entry in the data slot.  That makes the portfolio
                ## valid and will allow us to record missingness more
                ## accurately.

                curr.holdings <- expandData(curr.holdings)
                perf <- performance(curr.holdings, market.data = data)

                ## Update custom summary interface, if specified

                if(!is.null(object@summary.interface)){
                  sim.summary <- updateSummary(sim.summary,
                                               sim.data,
                                               sim.trades,
                                               curr.holdings,
                                               object@out.loc)
                }

                ## Create a result object.
               
                sp.result <- new("simResultSinglePeriod")

                ## Collect start instant data.
                
                sp.result@start.data@instant      <- periods$start[i]
                sp.result@start.data@equity.long  <-
                  portfolio:::mvLong(curr.holdings)
                sp.result@start.data@equity.short <-
                  portfolio:::mvShort(curr.holdings)
                sp.result@start.data@size.long  <-
                  portfolio:::sizeLong(curr.holdings)
                sp.result@start.data@size.short <-
                  portfolio:::sizeShort(curr.holdings)

                start.exposure <- exposure(curr.holdings, exp.var = object@exp.var)
                
                if(length(object@exp.var) > 0 && !is.null(start.exposure)){
                  sp.result@start.data@exposure <- start.exposure
                }
                if(length(object@contrib.var) > 0){
                  sp.result@period.data@contribution <-
                    contribution(curr.holdings,
                                 contrib.var = object@contrib.var,
                                 market.data = data)
                }
                
                
                sp.result@start.data@holdings <- curr.holdings


                trades <- new("trades", trades = trades)
                
                ## Collect period data.
                
                sp.result@period.data@period               <- period
                sp.result@period.data@turnover             <- turnover
                sp.result@period.data@universe.turnover    <- universe.turnover
                sp.result@period.data@performance          <- perf
                sp.result@period.data@trades               <- trades
                
                ## Expose holdings.
                
                curr.holdings <- expose(curr.holdings, trades)

                ## End instant data
                
                curr.holdings <-
                  updatePrices(curr.holdings, data$id, data$end.price)

                ## Calculate holdings weights

                curr.holdings <- calcWeights(curr.holdings)
                
                
                sp.result@end.data@instant      <- periods$end[i]
                sp.result@end.data@equity.long  <-
                  portfolio:::mvLong(curr.holdings)
                sp.result@end.data@equity.short <-
                  portfolio:::mvShort(curr.holdings)
                sp.result@end.data@size.long  <-
                  portfolio:::sizeLong(curr.holdings)
                sp.result@end.data@size.short <-
                  portfolio:::sizeShort(curr.holdings)

                if(length(object@exp.var) > 0){
                  sp.result@end.data@exposure <-
                    exposure(curr.holdings, exp.var = object@exp.var)
                }

                sp.result@end.data@holdings <- curr.holdings
                
                ## Saving the result object usually involves throwing
                ## away some information that we don't care about and
                ## would only slow us down.

                ## Here we save the pared down result and collect in
                ## our master result object.  The fully blown
                ## 'sp.result' object, however, is available to the
                ## next iteration if needed.  Should we be passing in
                ## to getSimTrades?
                
                sp.result.saved <- saveOut(sp.result,
                                           type    = object@out.type,
                                           out.loc = object@out.loc,
                                           verbose = verbose)

                ## The ith period's results get stored in the master
                ## result's data slot.  Preallocation prevents being
                ## out of bounds.
                
                master.result@data[[which(periods$period == period)]] <-
                  sp.result.saved

              }, silent = TRUE)

              ## If an error is encountered during the current
              ## iteration, make sure we revert to the previous set of
              ## holdings.  Perhaps a restart would be more elegant
              ## here.
              
              if(inherits(try.result, "try-error")){

                curr.holdings <- prev.holdings

                master.result@errors[[format(period)]] <-
                  as.character(try.result)
                
                ## User interrupt
                
                if(try.result[1] == ""){
                  stop("Simulation stopped by user.")
                }

                if(verbose){
                  cat("Error encountered during processing of period ",
                      format(period), ": ", try.result, "\n", sep = "")
                }
              }

            }

            master.result@summary.interface <- sim.summary

            saveOut(master.result, out.loc = object@out.loc, verbose = verbose)
            
            invisible(master.result)
          }
          )
