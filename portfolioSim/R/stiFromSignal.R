################################################################################
##
## $Id: stiFromSignal.R 1229 2007-10-02 14:07:55Z enos $
##
## Generate trades based on signal data stored in a data frame.  The
## way we do this is by determining a notion of target portfolio each
## period.  During the initial period we create the target using
## latest information as of that period.  In subsequent periods, if
## signal data is found for the current period a rebalance occurs.
## That is, we form a new target portfolio and use it as a basis for
## determining which trades to make.
##
################################################################################

setMethod("initialize",
          signature(.Object = "stiFromSignal"),
          function(.Object, ...){
            .Object <- callNextMethod()
            
            .Object@target <- new.env()
            .Object
          }
          )

setMethod("getSimTrades",
          signature(object         = "stiFromSignal",
                    period         = "orderable",
                    holdings       = "portfolio",
                    sim.data       = "simData",
                    verbose        = "logical"),
          function(object, period, holdings, sim.data, verbose){

            stopifnot(length(period) == 1)

            ## Create a new, empty input object

            sim.trades <- new("simTrades")

            ## Determine if we're rebalancing.  The client must supply
            ## periods on which to rebalance in the rebal.on slot.  We
            ## also implicitly "rebalance" if there is no target as a
            ## side effect of the period manipulation above.

            rebalance <- period %in% object@rebal.on
            data      <- sim.data@data

            ## Create holdings if there are none.
            
            if(is.null(holdings)){
              holdings <- new("portfolio", data = data)
            }

            ## Prices need to be set to the start-of-period prices so
            ## that weights and initial position market values are
            ## correct.  We are guaranteed this column in sim.data.

            target <- mget("target", envir = object@target,
                           inherits = FALSE , ifnotfound = list(NULL))$target
            
            if(rebalance || is.null(target)){

              stopifnot(object@in.var %in% names(data))

              ## Here's some hackery.  We should nail down how to
              ## handle NA prices when forming a portfolio.  Shares
              ## should never be NA!

              is.na(data[[object@in.var]]) <- is.na(data$start.price)

              target <- new("portfolio",
                            data      = data,
                            in.var    = object@in.var,
                            price.var = "start.price",
                            type      = object@type,
                            size      = object@size,
                            sides     = object@sides,
                            equity    = object@equity)

              # browser()
            }
            else{
              target@data <- data
              target      <- calcWeights(target)
            }

            assign("target", target, envir = object@target)
            
            if(isTRUE(all.equal(object@trading.style, "immediate"))){

              ## The holdings portfolio is already priced correctly by
              ## the sim.
              
              tl <- new("tradelist",
                        data      = data,
                        orig      = holdings,
                        target    = target,
                        type      = "all",
                        price.var = "start.price",
                        unrestricted = TRUE)

              ## Using style 'immediate' means returning orders to
              ## move completely from holdings to target.  Note the
              ## use of the unrestricted flag, providing unfettered
              ## trading.

              sim.trades@trades <- tl@final
            }
            else if(isTRUE(all.equal(object@trading.style, "percent.volume"))){

              stopifnot("volume" %in% names(data))

              tl <- new("tradelist",
                        data      = data,
                        orig      = holdings,
                        target    = target,
                        type      = "all",
                        price.var = "start.price",
                        unrestricted = TRUE
                        )

              x <- tl@final@trades

              ## Don't worry about side changes here.

              x <- x[order(x$id, match(x$side, c("C","S","B","X"))),]
              x <- x[!duplicated(x$id),]
              
              pct <- 15

              x$volume <- data$volume[match(x$id, data$id)]
              x$new.shares <- ifelse(x$shares > x$volume * (pct / 100),
                                     x$volume * (pct / 100),
                                     x$shares)
              
              x$shares <- x$new.shares
              x <- subset(x, shares > 0)
              
              sim.trades@trades <- new("trades", trades = x[c("id","side","shares")])

            }
            else if(isTRUE(all.equal(object@trading.style, "15.pct.vol.to.equity"))){

              data$round.lot <- 1
              data$a.6.s <- ifelse(is.na(data$alpha.6), 0, data$alpha.6)

              data$volume <- data$md.volume.120.d
              
              
              tl <- new("tradelist",
                        turnover  = object@turnover,
                        data      = data,
                        orig      = holdings,
                        target    = target,
                        type      = "ranks",
                        price.var = "start.price",
                        to.equity = TRUE,
                        tca       = c("volume.15.pct.only"),
                        rank.gain.min = -5000,
                        verbose = verbose,
                        chunk.usd = object@chunk.usd,
                        sorts     = list(a.6.s = 1)
                        )
              
              sim.trades@trades <- tl@final

            }  
            else{
              stop(paste("Unknown trading style", object@trading.style))
            }

            sim.trades@period <- period

            ## Sort for predictable output.

            if(nrow(sim.trades@trades@trades) > 0){
              sim.trades@trades@trades <-
                sim.trades@trades@trades[order(sim.trades@trades@trades$id),]
            }
            
            invisible(sim.trades)
          }
          )
