################################################################################
##
## $Id: tradelist.R 1311 2008-10-31 17:38:03Z enos $
##
## First pass at a trading system that enables daily trading
## at a given turnover rate.
##
## Constructs the tradelist with 7 methods. Stores the result of each
## step in a slot with a name similar to the method's name.
##
## 1. candidates - Determine set of candidate trades by "diffing"
##      original and target portfolios.  Removes candidates with external
##      restrictions.
##
## 2. ranks - Defines sorts on subsets of trades and order trades by
##      attractiveness within each sort.
##
## 3. chunks - Breaks each candidate into chunks of size
##      object@chunk.usd.  Computes trade-cost-adjusted rank for each
##      chunk.
##
## 4. swaps - Pair up B/S, C/X according to tca.rank quality.
##      "Dummy chunks" are created before the pairing process begins to
##      force the proper number of trades necessary to reach the target
##      equity on each side.
##
## 5. swaps.actual - Select the maximum number of swaps that are
##      allowed by object@turnover.  Order of selection is based on
##      total swap rank gain.
##
## 6. chunks.actual - Break each swap back into separate chunks,
##      removing dummy chunks.
##
## 7. actual - Aggregate chunks.
##
##
################################################################################

setMethod("initialize",
          signature(.Object = "tradelist"),
          function(.Object, orig, target, ...){

            .Object <- callNextMethod(.Object, ...)

            ## Allows creation of an "empty" tradelist
                                     
            if(!missing(orig) && !missing(target)){
              
              ## Target equity has no default in the prototype.  If
              ## none was supplied by the client in the call to 'new',
              ## use the equity of the portfolio 'target'.  If
              ## 'target' has no equity, use the average of the market
              ## values for the long and short sides.

              if(length(.Object@target.equity) == 0){
                if(length(target@equity) == 1){
                  .Object@target.equity <- target@equity
                }
                else{
                  mv.long <- mvLong(target)
                  mv.short <- mvShort(target)

                  if(mv.long > 0 && mv.short == 0){
                    .Object@target.equity <- mv.long
                  }
                  else if(mv.long == 0 && mv.short > 0){
                    .Object@target.equity <- abs(mv.short)
                  }
                  else{
                    .Object@target.equity <- mean(mv.long, abs(mv.short))
                  }
                }
              }

              .Object <- calcCandidates(.Object, orig, target)


              if(isTRUE(all.equal(.Object@type, "ranks"))){

                if(is.na(.Object@turnover) || .Object@turnover <= 0){
                  stop("Nothing to do because of turnover instruction")
                }
                if(isTRUE(.Object@unrestricted)){
                  stop("Cannot create tradelist of type 'ranks' with unrestricted set to TRUE")
                }

                
                if(is.null(.Object@sorts)){
                  warning("Using default.sort")
                }

                .Object <- calcRanks(.Object)
                .Object <- calcChunks(.Object)
                .Object <- calcSwaps(.Object)
                .Object <- calcSwapsActual(.Object)
                .Object <- calcChunksActual(.Object)
                .Object <- calcActual(.Object)
                .Object <- calcFinal(.Object)
                
                ## Need to add calcFinal to this sequence.
              }
              else if(isTRUE(all.equal(.Object@type, "all"))){

                ## Do all candidate trades.  I need to modify
                ## calcActual and calcFinal to do the right thing for
                ## the "all" type.  Right now this manually sets the
                ## contents of the final stlot.

                final.trades <- new("trades", trades =  .Object@candidates[finalCols(.Object)])
                .Object@final <- final.trades
              }
            }            
            .Object
          }
          )

setMethod("calcCandidates",
          signature(object = "tradelist",
                    orig   = "portfolio",
                    target = "portfolio"),
          function(object, orig, target){

            if(object@verbose) cat("Calculating candidates...\n")
            
            validObject(object)
 
            ## Copies the user-specified id column, id.var, to a
            ## column more easily referenced in the code, id

            if(object@id.var != "id"){
              object@data$id <- object@data[[object@id.var]]
            }

            ## We need to add a "symbol.var" slot to the tradelist
            ## object.  Before that time, we check to see whether
            ## symbol is a column in object@data.  If it is, we use it
            ## in row.names; o'wise just use id.

            if(!"symbol" %in% names(object@data)){
              object@data$symbol <- object@data$id
            }
            
            ## Merge the original and target portfolios, keeping
            ## rows where no match is found (exiting/entering
            ## positions).  Map NA -> 0, then keep only rows that
            ## will result in trades after adjusting each diff to be
            ## a round lot multiple.

            orig.positions   <- orig@shares
            target.positions <- target@shares
            
            if(nrow(orig.positions) == 0 &&
               nrow(target.positions) == 0){
              stop("Both portfolios are empty")
            }

            ## So that the merges work, trick the empty portfolio
            ## into thinking it has a single position (that also
            ## occurs in the other portfolio).  Since shares will be
            ## set to 0 we'll get the right trade for the copied
            ## candidate.
            
            if(nrow(orig.positions) == 0){
              orig.positions   <- target.positions[1,]
              orig.positions$shares <- 0
            }
            if(nrow(target.positions) == 0){
              target.positions <- orig.positions[1,]
              target.positions$shares <- 0
            }
            
            p.diff <- merge(orig.positions, target.positions,
                            by = "id", all = TRUE,
                            suffixes = c(".orig", ".target"))

            p.diff <- merge(p.diff, object@data,
                            by = "id", all.x = TRUE)

            p.diff$orig   <- ifelse(is.na(p.diff$shares.orig),
                                    0, p.diff$shares.orig)
            p.diff$target <- ifelse(is.na(p.diff$shares.target),
                                    0, p.diff$shares.target)

            ## Save off some important information about the original
            ## and target portfolios.

            object@mv.long.orig  <- mvLong(orig)
            object@mv.short.orig <- mvShort(orig)

            if(nrow(p.diff) > 0 && !"round.lot" %in% names(p.diff)){
              p.diff$round.lot <- 1
            }
           
            stopifnot(all(!is.na(p.diff$orig)))
            stopifnot(all(!is.na(p.diff$target)))

            ## We need to make application of round lots settable by a
            ## slot in the tradelist object.
            ##
            ## Each diff is forced into a round lot in this mapply.
            ## Closing trades are a special case and are not
            ## adjusted to the nearest round lot so that positions
            ## can be properly exited.
            
            p.diff$target <-
              mapply(
                     function(orig, target, round.lot){
                       if(target == 0){
                         return(target)
                       }

                       ## If there is no round lot, treat as 1.
                       
                       round.lot[is.na(round.lot)] <- 1
                       
                       current.diff <- target - orig
                       round.lot.diff <- .nearest.multiple(current.diff, round.lot)
                       orig + round.lot.diff
                     },

                     p.diff$orig,
                     p.diff$target,
                     p.diff$round.lot
                     )

            ## removes stocks from the tradelist that have the same
            ## number of shares in the original and target portfolios

            p.diff <- subset(p.diff, orig != target)

            ## Side changes are broken up here replacing each side change row
            ## with two rows, one for entry and one for exit.
            ##
            ## When sides are then calculated, we can be sure that each row
            ## corresponds to exactly one side in {"B","S","X","C"}.

            side.changes <- subset(p.diff, orig * target < 0)
            p.diff <- p.diff[! p.diff$id %in% side.changes$id,]

            if (nrow(side.changes) > 0) {

              exit  <- side.changes
              exit$target <- 0
              p.diff <- rbind(p.diff, exit)
              
              enter <- side.changes
              enter$orig <- 0
              p.diff <- rbind(p.diff, enter)
              
            }

            ## Calculates sides based on differences in original and target portfolios

            p.diff$side <-
              mapply(
                     function(orig, target) {
                       if(     orig >= 0 && target >= 0 && target > orig){
                         "B"
                       }
                       else if(orig >= 0 && target >= 0 && target < orig){
                         "S"
                       }
                       else if(orig <= 0 && target <= 0 && target > orig){
                         "C"
                       }
                       else if(orig <= 0 && target <= 0 && target < orig){
                         "X"
                       }
                       else{
                         stop("Encountered NO-OP, probably in error")
                       }
                     },
                     p.diff$orig, p.diff$target
                     )

            p.diff$shares <- abs(p.diff$target - p.diff$orig)
            p.diff$mv     <- (p.diff$target - p.diff$orig) * p.diff[[object@price.var]]
            
            if(!object@unrestricted){

              ## Basic data check.  Remove stocks that aren't in the
              ## data slot and haven't already been removed.  If
              ## trading is unrestriced, there is no use of the data
              ## slot so these steps are not performed.

              if(any(!p.diff$id %in% object@data$id)){
                object@restricted <-
                  rbind(object@restricted,
                        cbind(p.diff[!p.diff$id %in% object@data$id,],
                              data.frame(reason = "Not in data slot")
                              )[restrictedCols(object)])

                p.diff <- p.diff[p.diff$id %in% object@data$id,]
              }

              ## Strictly speaking, there is a set of data (columns)
              ## that are required for proper functioning of the
              ## tradelist calculation process.  Here we check for NA
              ## price and volume.  Having 0 for these values is OK.

              if(any(is.na(p.diff[[object@price.var]]))){
                object@restricted <-
                  rbind(object@restricted,
                        cbind(p.diff[is.na(p.diff[[object@price.var]]),],
                              data.frame(reason = "Missing price")
                              )[restrictedCols(object)])
                
                p.diff <- p.diff[!is.na(p.diff[[object@price.var]]),]
              }

              if(any(is.na(p.diff$volume))){
                object@restricted <-
                  rbind(object@restricted,
                        cbind(p.diff[is.na(p.diff$volume),],
                              data.frame(reason = "Missing volume")
                              )[restrictedCols(object)])
                
                p.diff <- p.diff[!is.na(p.diff$volume),]
              }
              
              ## Apply region restrictions.

              ## For now, only if the user supplies the proper column, region,
              ## and supplies a region upon which to filter in object@regions
              ## will any region filtering take place.  This feature is still
              ## under development.
              
              if(length(object@regions)){
                if("region" %in% names(object@data)){
                  p.diff <- p.diff[p.diff$region %in% object@regions,]
                }
                else{
                  warning("regions specified but no 'region' column found in slot 'data'")
                }
              }

              ## To simpify matters greatly, only include the closing portion
              ## of a side change.  This will make it easier to select trades
              ## from a given side to even out exposures.  The worst side
              ## effect is that it will take two days to implement the switch.
              ##
              ## After the trade is broken up in the iterative process,
              ## enter portions of the side change are pushed to the restricted
              ## data frame.

              if(nrow(p.diff) > 0){
                side.change.enter <- 
                  subset(p.diff,
                         id %in% subset(as.data.frame(table(p.diff$id,
                                                            dnn = "id")),
                                        Freq > 1)$id
                         & side %in% c("X","B"))

                p.diff <- 
                  subset(p.diff,
                         !(id %in% subset(as.data.frame(table(p.diff$id,
                                                              dnn = "id")),
                                          Freq > 1)$id
                           & side %in% c("X","B")))
                
                if(nrow(side.change.enter) > 0){
                  side.change.enter <- cbind(side.change.enter,
                                             data.frame(reason =
                                                        "Side change enter"))
                  object@restricted <-
                    rbind(object@restricted,
                          side.change.enter[restrictedCols(object)])
                }
              }
              
              ## Remove trades that violate any restrictions included in
              ## this object.  Save those restricted trades in the
              ## 'restricted' slot.

              if(nrow(object@restrictions) > 0){
                restricted <- merge(p.diff, object@restrictions,
                                    by.x = c("id","side"),
                                    by.y = c("id","type"))
                
                object@restricted <-
                  rbind(object@restricted,
                        restricted[restrictedCols(object)])
                p.diff <- subset(p.diff, ! id %in% restricted$id)
              }

              ## Enforce trade market value minimum.  Ideally what we
              ## want to do is enforce the trade market value minimum
              ## before we even break up side changes.  As an exception,
              ## we always want to perform enters and exits from the
              ## portfolio.

              if(nrow(p.diff[abs(p.diff$mv) < object@trade.usd.min &
                             p.diff$orig   != 0 &
                             p.diff$target != 0,]) > 0){

                object@restricted <-
                  rbind(object@restricted,
                        cbind(p.diff[abs(p.diff$mv) < object@trade.usd.min &
                                     p.diff$orig   != 0 &
                                     p.diff$target != 0,],
                              data.frame(reason = "Trade too small")
                              )[restrictedCols(object)])

                p.diff <- p.diff[!(abs(p.diff$mv) < object@trade.usd.min &
                                   p.diff$orig   != 0 &
                                   p.diff$target != 0),]
              }
              
              
              ## After every merge (and perhaps after every rbind, I
              ## have to check) we lose row names.  Here I add them
              ## back to the list of restricted items.  If we were to
              ## add to this list in other functions, I presume we'd
              ## have to set them again.

              rowname.tag <- object@data$symbol[match(object@restricted$id,
                                                      object@data$id)]
              rowname.tag <- ifelse(is.na(rowname.tag), object@restricted$id,
                                    rowname.tag)
                      
              row.names(object@restricted) <-
                paste(rowname.tag,
                      object@restricted$side,
                      object@restricted$reason, sep = " ")
              
            }

            ## If we're running in unrestricted mode, there can be
            ## more than one row per security, as in a side change.

            ## Do we require that symbol is present?  Should we have a
            ## symbol.var?
            
            p.diff$symbol <-
              ifelse(is.na(p.diff$symbol), p.diff$id, p.diff$symbol)
            
            if(any(duplicated(p.diff$symbol))){
              if(any(duplicated(paste(p.diff$symbol, p.diff$side)))){
                row.names(p.diff) <- paste(p.diff$id, p.diff$side)
              }
              else{
                row.names(p.diff) <- paste(p.diff$symbol, p.diff$side)
              }
            }
            else{
              row.names(p.diff) <- p.diff$symbol
            }

            p.diff            <- p.diff[order(row.names(p.diff)),]
            object@candidates <- p.diff[candidatesCols(object)]

            object

          }
          )

setMethod("calcRanks",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating ranks...\n")
            
            validObject(object)

            if(nrow(object@candidates) == 0){
              warning("No candidates found, doing nothing.")
              return(object)
            }

            ## Makes sure that the candidates data frame has been
            ## correctly formed and contains the necessary columns

            req.cols <- c("id", "side", "shares", "mv")
            
            if(!all(req.cols %in% names(object@candidates))){
              stop("\n\"candidates\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@candidates)], sep = "")
            }

            ## Makes sure that there is an "id" column in
            ## object@data.

            stopifnot("id" %in% names(object@data))

            ## Checks that the user has supplied in "data" all the
            ## sorts specified in the sorts list.

            if(!all(names(object@sorts) %in% names(object@data))){
              stop("Not all sorts have a column in the 'data' slot")
            }
            
            ranks <- merge(object@candidates, object@data, by = "id")

            stopifnot(nrow(ranks) > 0)

            ## This should be replaced with ranks$object@symbol once
            ## we add a symbol.var slot
            
            row.names(ranks) <- ranks$symbol
            
            ## Sort a data frame of "trades", where a trade is a side
            ## token in {B,S,C,X}, by a numeric value.  More positive
            ## values are better for B,C, and more negative are better
            ## for S,X.  First we rank within side. "Interleaves" each
            ## side's ranked orders using the following order: B, S,
            ## C, X

            .sort.trades <- function(x, by.var, side.var = "side"){
              cols.orig <- names(x)

              ## (1) Splits a data frame into a list of data frames by
              ## side (B,S,C,X). (2) ranks the trades within each side
              ## such that each trade has a rank (1, 2, 3...nrow(x))
              ## (3) Combines the data frames into a single data frame
              ## using "unsplit"

              x$rank <-
                unsplit(lapply(split(x, x[[side.var]]),
                               function(x){
                                 rank(ifelse(x$side %in% c("B","C"), -1, 1) * x[[by.var]],
                                      na.last = TRUE, ties.method = "first")
                               }
                               ), x[[side.var]])
              
              ## Unsplit leaves us with up to 4 trades that share a
              ## single rank. Sorts the trades by rank s.t. when
              ## trades of different sides share a rank, we always
              ## order them as follows: B, S, C, X

              x <- x[order(x$rank, match(x$side, c("B","S","C","X"))),]
              x[cols.orig]
            }

            rank.sorts        <- list()
            sorts             <- object@sorts
            
            ## Applies user-defined sorts to generate ranks

            for(s in names(sorts)){
              sorted            <- .sort.trades(ranks, s)
              sorted            <- sorted[!is.na(sorted[[s]]),]

              ## If there are no orders in this sort, remove it from
              ## the list of sorts in this object.

              if(nrow(sorted) == 0){
                object@sorts <- object@sorts[names(object@sorts) != s]
                next()
              }

              ## Applies relative weights assigned to sorts. For
              ## example, default.sort = 1, but another sort, "foo",
              ## might be defined foo = 1/10

              sorted$rank <- 1:nrow(sorted) / sorts[[s]]

              ## sorted$symbol should be sorted$object@symbol
              
              row.names(sorted)      <- sorted$symbol
              rank.sorts[[s]]        <- sorted[c(candidatesCols(object), s, "rank")]
            }

            ## Save off the list of all rank sorts.
            
            object@rank.sorts <- rank.sorts

            ## Associate with each trade another rank.  For each stock
            ## we take the most attractive rank among all sorts.  This
            ## rank will not be unique in the case where a two trades
            ## on the same side have been given the same rank by
            ## different sorts.
            
            ## Stacks the sorts in a data frame, with the sorts with
            ## lower indices in the "sorts" slot appearing closer to
            ## the top of the data frame
            
            all.sorts  <- do.call(rbind, lapply(rank.sorts, function(x) { x[c("id","rank")] }))

            ## Removes duplicates from "all.sorts". Associates with
            ## each trade the best rank any sort has assigned it

            top.ranks  <- aggregate(all.sorts[c("rank")], by = list(id = all.sorts$id), min)

            ## Updates the trade rankings using the latest measure of
            ## rank as determined in the two prior steps

            ranks$rank <- top.ranks$rank[match(ranks$id, top.ranks$id)]
            
            ## Synthesise another rank, "rank.t" to use when ranking
            ## all of our trades with a number throughout the rest of
            ## the process.  "rank.t" will not be unique if there are
            ## the same number of trades within the buy and cover
            ## sides or the sells and shorts sides

            ## 'rank.ws' means 'rank within side'.  We're scaling
            ## uniformly spaced ranks on [0,1] to percentiles of a
            ## normal distribution.  Buys and covers, therefore, will
            ## have increasing rank as scaled rank approaches 1;
            ## shorts and sells have increasing rank as scaled rank
            ## approaches 0.  Covers and buys with worse ranks get
            ## "better" ranks because we multiply the ranks of covers
            ## and buys by -1.  Since lesser ranks are better, the
            ## buys and covers that have initially been ranked poorly
            ## (greater numbers like 100, 101, 102), now have better
            ## ranks (-100, -101, -102).

            ranks$rank.ws <-
              unsplit(lapply(split(ranks, ranks$side),
                             function(x) { rank(ifelse(x$side %in% c("C","B"), -1, 1) * x$rank,
                                                na.last = "keep") } ), ranks$side)

            ## Determines the maximum rank within each side and adds 1
            ## so we don't ever pass a value of 1 to "qnorm"

            r.max  <- tapply(ranks$rank.ws, ranks$side, max, na.rm = TRUE) + 1

            r.mult <- unlist(list(X = 0.15, C = 0.15, S = 0.15, B = 0.15))
            r.add  <- unlist(list(X = 0,    C = 0.85, S = 0,    B = 0.85))

            ## Scales buys and covers to values [0.85-1), and sells and shorts, (0, 0.15]

            ranks$rank.ws <-
              as.vector(r.mult[ranks$side] * ranks$rank.ws / r.max[ranks$side] +
                        r.add[ranks$side])
            
            ## Maps buys and covers to the normal distribution

            ranks$rank.t  <- qnorm(ranks$rank.ws)

            object@ranks  <- ranks[ranksCols(object)]
            object
          }
          )

setMethod("calcChunks",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating chunks...\n")
            
            validObject(object)

            if(nrow(object@ranks) == 0){
              warning("No ranks found, doing nothing.")
              return(object)
            }
            
            ## Makes sure that the ranks data frame has been
            ## correctly formed and contains the necessary columns

            req.cols <- c("id", "side", "shares", "mv", "rank.t", names(object@sorts))
            
            if(!all(req.cols %in% names(object@ranks))){
              stop("\n\"ranks\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@ranks)], sep = "")
            }

            stopifnot(c("id", "volume") %in% names(object@data))
            
            ## Also, removes the sort columns so that the merge does
            ## not contain two of the same columns with different
            ## suffixes.
            
            candidates <- merge(object@ranks,
                                object@data[!names(object@data) %in% names(object@sorts)], by = "id")

            ## The chunks slot contains a data frame with each trade
            ## broken up into reasonably-sized chunks.  Each chunk
            ## may have a different trade-cost-adjusted rank
            ## assuming that the first chunk is cheaper to do than
            ## the second, etc.
            ##
            ## This is an iterative process and involves making sure
            ## each chunk is a round lot multiple and that chunk
            ## remainders are handled properly.
            ##
            ## We preallocate space for storing chunks to avoid a
            ## bottleneck here.
            
            chunks <- NULL
            did.prealloc <- FALSE

            ## Calculate an upper bound on the number of chunks we'll
            ## have.

            stopifnot(all(!is.na(candidates$mv)))
            max.num.chunks <- sum(ceiling(abs(candidates$mv) / object@chunk.usd))
            chunks.pos <- 1

            ## Slim down candidates into candidates.slim.  It contains
            ## the minimum amount of information required to make
            ## chunks.

            stopifnot("round.lot" %in% names(candidates))
            
            slim.cols <- c("id", object@price.var, "round.lot", "shares")
            candidates.slim <- candidates[slim.cols]

            ## Everything else, which we'll merge in later, goes in
            ## candidates.supp.

            supp.cols <-
              c("id",names(candidates)[!names(candidates) %in% names(candidates.slim)])
            candidates.supp <- candidates[supp.cols]
            
            if(nrow(candidates.slim) > 0){
              for(i in 1:nrow(candidates.slim)){

                ## Each chunk starts out as a copy of its
                ## corresponding candidate row.  Chunk numbering
                ## begins at 1, and chunk size (in shares) must be
                ## calculated separately for each candidate.
                
                p.row <- candidates.slim[i,]
                chunk <- 1
                chunk.shares <- .nearest.multiple(object@chunk.usd / p.row[[object@price.var]],
                                                              p.row$round.lot)

                ## If chunk.shares is 0, that is, a single chunk isn't
                ## half a round lot and has been rounded to 0, set the
                ## chunk.shares to the total shares of the candidate.

                if(chunk.shares == 0){
                  chunk.shares <- p.row$shares
                }
                
                total.shares <-  p.row$shares

                while(total.shares > 0){

                  ## After copying the chunk, increment the chunk
                  ## number.  The portion of the candidate trade
                  ## that falls under the chunk amount becomes its
                  ## own chunk.
                  
                  p.chunk              <- p.row
                  p.chunk$chunk.shares <- min(chunk.shares, total.shares)
                  total.shares         <- total.shares - p.chunk$chunk.shares
                  p.chunk$chunk        <- chunk
                  chunk                <- chunk + 1

                  if(!did.prealloc){
                    chunks <- p.chunk[FALSE,]
                    chunks <- chunks[1:max.num.chunks,]
                    did.prealloc <- TRUE
                  }

                  chunks[chunks.pos,] <- p.chunk
                  chunks.pos <- chunks.pos + 1

                }
              }
            }

            if(chunks.pos <= nrow(chunks)){
              chunks <- chunks[1:(chunks.pos - 1),]
            }
            
            ## Now put back all the information we took out.

            chunks <- merge(chunks, candidates.supp, by = "id", all.x = TRUE)
            
            
            ## Record the (signed) market value of each chunk, and
            ## calculate our measure of trading volume for the
            ## security corresponding to each chunk.
            
            chunks$chunk.mv <- chunks$chunk.shares *
              chunks[[object@price.var]] * sign(chunks$mv)
            
            chunks$trading.volume <- chunks$volume

            ## The cumulative sum of chunk market value is set here.
            ## This value is used to determine the percentile range
            ## of the day's volume this chunk represents, which in
            ## part determines our measure of trade-cost-adjusted
            ## rank.
            
            chunks <- chunks[order(chunks$id, chunks$chunk),]
            chunks$cumsum.chunk.shares <- 
              unsplit(lapply(split(chunks$chunk.shares, chunks$id), cumsum),
                      chunks$id)

            stopifnot(all(!is.na(chunks$trading.volume)))

            ## Trade-cost-adjusted rank requires the five
            ## parameters listed below from each chunk.  See the
            ## .tca function in the utility function section below
            ## for more information.

            ## This is really ugly.  Need to refactor.
            
            chunks$tca.rank <- chunks$rank.t
            
            if("volume" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.volume,
                                         chunks$tca.rank,
                                         chunks$chunk.shares,
                                         chunks$cumsum.chunk.shares,
                                         chunks$side,
                                         chunks$trading.volume,
                                         chunks$chunk)
            }
            else if("volume.conservative" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.volume.conservative,
                                         chunks$tca.rank,
                                         chunks$chunk.shares,
                                         chunks$cumsum.chunk.shares,
                                         chunks$side,
                                         chunks$trading.volume,
                                         chunks$chunk)
            }
            else if("volume.15.pct.only" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.volume.15.pct.only,
                                         chunks$tca.rank,
                                         chunks$chunk.shares,
                                         chunks$cumsum.chunk.shares,
                                         chunks$side,
                                         chunks$trading.volume,
                                         chunks$chunk)
            }
            else if("volume.10.pct.only" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.volume.10.pct.only,
                                         chunks$tca.rank,
                                         chunks$chunk.shares,
                                         chunks$cumsum.chunk.shares,
                                         chunks$side,
                                         chunks$trading.volume,
                                         chunks$chunk)
            }
            if("entry" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.entry,
                                         chunks$tca.rank,
                                         chunks$side,
                                         chunks$orig,
                                         chunks$chunk)
            }
            if("too.short" %in% object@tca){
              chunks$tca.rank <- mapply(.tca.too.short,
                                         chunks$tca.rank,
                                         chunks$side,
                                         chunks$mv,
                                         object@mv.long.orig)
            }

            
            ## Row names for chunk are symbol + chunk number.  Chunk
            ## number 2 for IBM would be "IBM.2".
            
            row.names(chunks) <- paste(chunks$symbol, chunks$chunk, sep = ".")
            object@chunks <- chunks[chunksCols(object)]

            object
          }
          )


setMethod("calcSwaps",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating swaps...\n")
            
            validObject(object)

            if(nrow(object@chunks) == 0){
              warning("No chunks found, doing nothing.")
              return(object)
            }
            
            ## Makes sure that the chunks data frame has been
            ## correctly formed and contains the necessary columns

            req.cols <- c("id", "side", "shares", "mv", "rank.t",
                          "tca.rank", "chunk.shares", "chunk.mv",
                          "chunk", names(object@sorts))
            
            if(!all(req.cols %in% names(object@chunks))){
              stop("\n\"chunks\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@chunks)], sep = "")
            }

            ## Begin by splitting chunks into separate data frames by
            ## side.  Logic for dealing with long/short sides is
            ## specific to each side, so we don't mind this loss of
            ## generality.  Also, ensure that there are no null side
            ## buckets (only empty data frames).

            chunks <- object@chunks
            chunks$dummy.quality <- NA

            ## At this point, remove all chunks that have an NA
            ## tca.rank.  By this mechanism we can categorically
            ## exclude some chunks in our tca functions.

            chunks <- subset(chunks, !is.na(tca.rank))

            ## Split up chunks by side.
            
            split.chunks <- split(chunks, chunks$side)

            ## If any one of the types of trades are missing, set the
            ## data fame for that type (buy, sell, etc..) to have 0
            ## rows and the same columns as the data frame in the
            ## chunks slot

            if(is.null(split.chunks$B))
              buys <- chunks[FALSE,]
            else buys <- split.chunks$B
            
            if(is.null(split.chunks$S))
              sells <- chunks[FALSE,]
            else sells <- split.chunks$S

            if(is.null(split.chunks$C))
              covers <- chunks[FALSE,]
            else covers <- split.chunks$C

            if(is.null(split.chunks$X))
              shorts <- chunks[FALSE,]
            else shorts <- split.chunks$X
            
            to.equity.long  <-   object@target.equity - object@mv.long.orig
            to.equity.short <- - object@target.equity - object@mv.short.orig


            if(object@to.equity) {

              ## Creation of dummy chunks in order to get to target
              ## equity.

              ## Dummy chunks force the first traded dollars to
              ## address over/under exposure in respect to our target
              ## equity.  For example, if we are over invested by
              ## $20,000 on the long side, we introduce $20,000 worth
              ## of dummy buy chunks that look exceptionally
              ## attractive.  We won't trade these dummy chunks, but
              ## we will trade their paired sells. 

              ## Long side: over  = create sells
              ##            under = create buys

              if(to.equity.long > 0){
                num.dummy <-
                  sum(cumsum(buys$chunk.mv[order(buys$tca.rank,
                                                 decreasing = TRUE)])
                      < to.equity.long)
                dummy.sells <- dummyChunks(object, "S", num.dummy, "good")
                sells <- rbind(sells, dummy.sells)
              }
              else if(to.equity.long < 0){
                num.dummy <-
                  sum(cumsum(sells$chunk.mv[order(sells$tca.rank,
                                                  decreasing = FALSE)])
                      > to.equity.long)
                dummy.buys  <- dummyChunks(object, "B", num.dummy, "good")
                buys <- rbind(buys, dummy.buys)
              }

              ## Short side: over  = create shorts
              ##             under = create covers

              if(to.equity.short > 0){
                num.dummy <-
                  sum(cumsum(covers$chunk.mv[order(covers$tca.rank,
                                                   decreasing = TRUE)])
                      < to.equity.short)
                dummy.shorts <- dummyChunks(object, "X", num.dummy, "good")
                shorts <- rbind(shorts, dummy.shorts)
              }
              else if(to.equity.short < 0){
                num.dummy <-
                  sum(cumsum(shorts$chunk.mv[order(shorts$tca.rank,
                                                   decreasing = FALSE)])
                      > to.equity.short)
                dummy.covers <- dummyChunks(object, "C", num.dummy, "good")
                covers <- rbind(covers, dummy.covers)
              }
            }

            ## Now we need to even up the number of sells and
            ## covers, so that trades don't get lost by the wayside.
            ## We do this by adding not very attractive trades that
            ## are sure to sink to the bottom of the list.

            if(nrow(buys) > nrow(sells))
              sells <-   rbind(sells,
                               dummyChunks(object, "S",
                                            nrow(buys) - nrow(sells),
                                            "bad"))
            
            if(nrow(buys) < nrow(sells))
              buys <-   rbind(buys,
                              dummyChunks(object, "B",
                                           nrow(sells) - nrow(buys),
                                           "bad"))
            if(nrow(covers) > nrow(shorts))
              shorts <- rbind(shorts,
                              dummyChunks(object, "X",
                                           nrow(covers) - nrow(shorts),
                                           "bad"))

            if(nrow(covers) < nrow(shorts))
              covers <- rbind(covers,
                              dummyChunks(object, "C",
                                           nrow(shorts) - nrow(covers),
                                           "bad"))

            ## Creating swaps.

            ## Swaps are formed by merging sorted buys/sells and
            ## covers/shorts.  The most attractive buys, for
            ## example, should be paired with the most attractive
            ## sells (that is, sells where tca.rank is lowest).
            ## There should be no more swaps than the minimum number
            ## of trades in each of the two sides involved.

            ## To achieve this effect, I sort the trades by
            ## tca.rank according to side, and then set row.names
            ## to be the row number.  The (implicit) all = FALSE of
            ## the merge makes sure we get the correct number of
            ## swaps.
            
            buys <- buys[order(buys$tca.rank, decreasing = TRUE),]
            if(nrow(buys) > 0){
              row.names(buys) <- 1:nrow(buys)
            }

            sells <- sells[order(sells$tca.rank, decreasing = FALSE),]
            if(nrow(sells) > 0){
              row.names(sells) <- 1:nrow(sells)
            }

            ## The 'enter' and 'exit' suffixes allow us to rbind all
            ## swaps at a later stage.

            swaps.long <- NULL
            if(nrow(buys) > 0 && nrow(sells) > 0){
              swaps.long <- merge(buys, sells, by = "row.names",
                                  suffixes = c(".enter",".exit"))

              ## The rank.gain of a swap on the long side is the
              ## simple swaperence of tca.rank.
            
              swaps.long$rank.gain <- (swaps.long$tca.rank.enter -
                                        swaps.long$tca.rank.exit)
            }
            
            covers <- covers[order(covers$tca.rank, decreasing = TRUE),]
            if(nrow(covers) > 0){
              row.names(covers) <- 1:nrow(covers)
            }

            shorts <- shorts[order(shorts$tca.rank, decreasing = FALSE),]
            if(nrow(shorts) > 0){
              row.names(shorts) <- 1:nrow(shorts)
            }

            swaps.short <- NULL
            if(nrow(shorts) > 0 && nrow(covers) > 0){
              swaps.short <- merge(shorts, covers, by = "row.names",
                                   suffixes = c(".enter",".exit"))

              ## The rank.gain measure on the short side requies the
              ## opposite sign from the long side.
            
              swaps.short$rank.gain <- - (swaps.short$tca.rank.enter -
                                           swaps.short$tca.rank.exit)
            }
            
            swaps <- rbind(swaps.long, swaps.short)
            swaps <- swaps[order(swaps$rank.gain, decreasing = TRUE),]

            ## It is not uncommon, particularly when dealing with
            ## small portfolios, to have dummy chunks lined up with
            ## each other at this point.  Remove them now.

            swaps <- swaps[!(swaps$id.enter == .dummy.id() &
                             swaps$id.exit  == .dummy.id()),]

            ## Row names for swaps are the chunk names separated by
            ## a column.  Since I don't keep around unnecessary
            ## information like symbol for easy merging with secref
            ## and other supplementary data sources, we have to
            ## recreate the chunk names from scratch.

            row.names(swaps) <-
              paste(paste(object@data$symbol[match(swaps$id.enter,
                                                     object@data$id)],
                          swaps$chunk.enter, sep = "."),
                    paste(object@data$symbol[match(swaps$id.exit,
                                                     object@data$id)],
                          swaps$chunk.exit, sep = "."),
                    sep = ",")
            
            object@swaps <- swaps
            
            object
            
          }
          )

setMethod("calcSwapsActual",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating actual swaps...\n")

            validObject(object)

            if(nrow(object@swaps) == 0){
              warning("No swaps found, doing nothing.")
              return(object)
            }
            
            ## Makes sure that the swaps data frame has been correctly
            ## formed and contains the necessary columns

            req.cols <- c("chunk.mv.enter", "chunk.mv.exit")
            
            if(!all(req.cols %in% names(object@swaps))){
              stop("\n\"swaps\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@swaps)], sep = "")
            }

            swaps <- object@swaps

            ## Remove all swaps that don't have a positive enough rank
            ## gain.  Note the strict inequality here, meaning if rank
            ## gain is -Inf the swap won't be done even if
            ## rank.gain.min == -Inf.

            swaps <- swaps[swaps$rank.gain > object@rank.gain.min,]

            ## Controlling turnover is actually the easy part.
            
            ## We're allowed to select swaps to trade until the
            ## unsigned dollar amount of the trades in the swaps
            ## exceeds our target turnover.  Note that it is
            ## important that dummy chunks have a zero market value
            ## so that they play nicely with this process.
            
            swaps.actual <- subset(swaps, cumsum(abs(chunk.mv.enter) +
                                                 abs(chunk.mv.exit)) <= object@turnover)

            ## If to.equity == TRUE, however, we don't want to worsen
            ## exposure by doing the last set of remaining trades.
            ## For each side we first want to do swaps with dummy
            ## chunks intended to bring us closer to target.equity,
            ## then any swaps not including dummy chunks.

            ## When we computed swaps, if to.equity == TRUE, we tagged
            ## each swap with a T/F flag indicating whether it did, in
            ## fact, help us get to the target equity as described
            ## above.  As far as our representation, this means we
            ## don't want to do any swap that involves a chunk with
            ## dummy.quality %in% "bad".

            if(object@to.equity){
              swaps.actual <- subset(swaps.actual, ! (dummy.quality.enter %in% "bad" | dummy.quality.exit %in% "bad"))
            }
            
            object@swaps.actual <- swaps.actual

            object
          }
          )

setMethod("calcChunksActual",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating actual chunks...\n")            

            validObject(object)

            if(nrow(object@swaps.actual) == 0){
              warning("No swaps.actual found, doing nothing.")
              return(object)
            }
            
            ## Makes sure that the swaps.actual data frame has been
            ## correctly formed and contains the necessary columns

            req.cols <- c("chunk.mv.enter", "chunk.mv.exit")

            if(!all(req.cols %in% names(object@swaps.actual))){
              stop("\n\"swaps.actual\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@swaps.actual)], sep = "")
            }

            swaps.actual <- object@swaps.actual

            ## Now that we've selected which swaps we're actually
            ## going to trade, we need to turn those swaps back into
            ## chunks.  Since enter/exit chunks within swaps can be
            ## identified by their column suffixes ('enter' and
            ## 'exit'), I use grep and sub to extract chunks.
            
            cols.enter <- grep("\\.enter", names(swaps.actual), value = TRUE)
            chunks.enter <- swaps.actual[cols.enter]
            names(chunks.enter) <- sub("\\.enter", "", names(chunks.enter))
            chunks.enter <- subset(chunks.enter, id != .dummy.id())

            cols.exit <- grep("\\.exit", names(swaps.actual), value = TRUE)
            chunks.exit <- swaps.actual[cols.exit]
            names(chunks.exit) <- sub("\\.exit", "", names(chunks.exit))
            chunks.exit <- subset(chunks.exit, id != .dummy.id())

            chunks.actual <- rbind(chunks.enter, chunks.exit)

            row.names(chunks.actual) <-
              paste(object@data$symbol[match(chunks.actual$id,
                                               object@data$id)],
                    chunks.actual$chunk, sep = ".")
            
            object@chunks.actual <- chunks.actual

            object

          }
          )

setMethod("calcActual",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating actual...\n")

            validObject(object)

            if(nrow(object@chunks.actual) == 0){
              warning("No chunks.actual found, doing nothing.")
              return(object)
            }
            
            ## Makes sure that the chunks.actual data frame has been
            ## correctly formed and contains the necessary columns

            req.cols <- c("id", "side", "shares", "mv", "rank.t",
                          "tca.rank", "chunk.shares", "chunk.mv",
                          "chunk", names(object@sorts))

            if(!all(req.cols %in% names(object@chunks.actual))){
              stop("\n\"chunks.actual\" data frame requires column(s) for ",
                   req.cols[!req.cols %in% names(object@chunks.actual)], sep = "")
            }

            ## Once we have a list of chunks we're ready to trade,
            ## we need to roll these chunks up on a per-security
            ## basis.
            
            actual <- subset(object@ranks,
                             id %in% object@chunks.actual$id)

            ## In future versions, this should be abstracted to "by =
            ## id.var".  Also, the strange subsetting operation that
            ## removes duplicate rows with the same names as the
            ## objects in object@sorts
            
            actual <- merge(actual, object@data[,!names(object@data)
                                                %in% names(object@sorts)], by = "id")

            actual$shares <- 0

            ## Using a tapply, sum up shares across securities.
            ## Perhaps a better way of linking this data is using a
            ## 'match' instead of creating a data frame, but I'm not
            ## sure how this would work when dealing with a named
            ## list instead of a data frame.
            
            shares.sum <-
              tapply(object@chunks.actual$chunk.shares,
                     object@chunks.actual$id, sum)
            shares.sum <- data.frame(id = names(shares.sum),
                                     sum = as.numeric(shares.sum))
            actual <- merge(actual, shares.sum, by = "id")
            actual$shares <- actual$sum

            actual$mv <-
              mapply(
                     function(side,shares,price){
                       mv <- shares * price
                       if (side %in% c("X","S"))
                         return(-1 * mv)
                       else
                         return(mv)
                     },
                     actual$side,
                     actual$shares,
                     actual[[object@price.var]]
                     )
            
            row.names(actual) <- actual$symbol
            actual <- actual[order(row.names(actual)),]
            
            object@actual <- actual[actualCols(object)]
            object
            
          }
          )

setMethod("calcFinal",
          "tradelist",
          function(object){

            if(object@verbose) cat("Calculating final...\n")

            validObject(object)

            if(nrow(object@actual) == 0){
              warning("No orders found in actual slot. Doing nothing.")
              return(object)
            }

            final.trades <- new("trades", trades = object@actual[finalCols(object)])
            
            object@final <- final.trades
            object
            
          }
          )

## The following methods are used to select the important columns
## from intermediate data frames for storage in the slots of this
## object.

setMethod("candidatesCols",
          "tradelist",
          function(object){
            c("id","orig","target","side","shares","mv")
          }
          )

setMethod("ranksCols",
          "tradelist",
          function(object){
            c(candidatesCols(object), names(object@sorts), "rank.t")
          }
          )

setMethod("actualCols",
          "tradelist",
          function(object){
            c("id","side","shares","mv", names(object@sorts), "rank.t")
          }
          )

setMethod("finalCols",
          "tradelist",
          function(object){
            c("id","side","shares")
          }
          )

setMethod("chunksCols",
          "tradelist",
          function(object){
            c(ranksCols(object), "tca.rank", "chunk.shares", "chunk.mv", "chunk")
          }
          )

setMethod("restrictedCols",
          "tradelist",
          function(object){
            c(candidatesCols(object), "reason")
          }
          )

setMethod("trimSide",
          "tradelist",
          function(object, side, value){
            trades.side <- object@actual[object@actual$side == side,]
            decreasing  <- side %in% c("S","X")
            trades.side <- trades.side[order(trades.side$rank.t, decreasing = decreasing),]
            trades.side <- trades.side[cumsum(abs(trades.side$mv)) < abs(value),]

            object@actual <- subset(object@actual, ! id %in% trades.side$id)
            object
          }
          )

## The show method displays the important slot information of this
## object in addition to views on the final trades from several
## perspectives.

setMethod("show",
          signature(object = "tradelist"),
          function(object) {

            if(length(object@turnover) == 0){
              object@turnover <- 0
              warning("turnover missing; setting to 0 for display")
            }
            if(length(object@target.equity) == 0){
              object@target.equity <- 0
              warning("target.equity missing; setting to 0 for display")
            }
            
            cat("Tradelist object\n\nInput:\n\n",
                sprintf("%-15s %12s", "Turnover:",      object@turnover),"\n",
                sprintf("%-15s %12s", "Target Equity:", object@target.equity), "\n",
                "\n\n",
                sep = ""
                )

            if(nrow(object@candidates) > 0){
              
              ## Input summary: total dollars to do, including
              ## restricted (all) and excluding restricted (doable).

              if(nrow(object@restricted) != 0){
                all.trades    <- rbind(object@candidates, object@restricted[candidatesCols(object)])
              }
              else{
                all.trades    <- object@candidates
              }
              
              doable.trades <- object@candidates

              ## Below, pct.close refers to how much trading we have
              ## left to do until we reach the target.  It only applies
              ## to the region considered in this tradelist; we assume
              ## that the book is 50/50 domestic/international to
              ## compute the denominator.
              
              cat(
                  sprintf("%-7s %6s %12s", "", "count", "value"), "\n",
                  sprintf("%-7s %6s %12s",
                          "all",
                          prettyNum(round(nrow(all.trades)), big.mark = ","),
                          prettyNum(round(sum(abs(all.trades$mv), na.rm = TRUE)), big.mark = ",")), "\n",
                  sprintf("%-7s %6s %12s",
                          "doable",
                          prettyNum(round(nrow(doable.trades)), big.mark = ","),
                          prettyNum(round(sum(abs(doable.trades$mv))), big.mark = ",")), "\n\n",
                  sep = "")

              ## Input summary: dollars required per day, using 20 day
              ## median volume only, limit of 10 days.
              
              doable.trades <- merge(doable.trades, object@data, by = "id", all.x = TRUE)
              doable.trades$days.required <- abs(doable.trades$shares) / (0.15 * doable.trades$volume)
              
              value.per.day <- function(x, day) {
                ifelse(x$days.required > day,
                       x$mv / x$days.required,
                       ifelse(day >= x$days.required & x$days.required > day-1,
                              (x$mv / x$days.required) * (x$days.required - floor(x$days.required)), NA)) }

              for (i in c(1:10)) doable.trades[paste("day.",i, sep = "")] <- value.per.day(doable.trades,i)

              by.day <- data.frame(value = sapply(doable.trades[grep('^day\\.\\d+$', names(doable.trades), perl = TRUE)],
                                     function(x) { sum(abs(x), na.rm = TRUE) }))

              mapply(function(x,y){ cat(sprintf("%-7s %6s %12s", x, "", prettyNum(round(y), big.mark = ",")), "\n") },
                     row.names(by.day),
                     by.day$value)
              
            }

            if(nrow(object@actual) > 0){
              
              sides <- split(object@actual, object@actual$side)

              cat("\nActual trades:\n\n")
              cat(
                  sprintf("%-7s %6s %12s", "side", "count", "value"), "\n",
                  sep = "")

              all.sides <- c("B","S","C","X")
              for(side in all.sides[all.sides %in% names(sides)]){
                s <- sides[[side]]
                cat(
                    sprintf("%-7s %6s %12s",
                            side,
                            prettyNum(round(nrow(s)), big.mark = ","),
                            prettyNum(round(sum(s$mv)), big.mark = ",")), "\n",
                    sep = "")
              }

              long  <- subset(object@actual, side %in% c("B","S"))
              short <- subset(object@actual, side %in% c("C","X"))
              cat(
                  sprintf("%-7s %6s %12s",
                          "\nnet (l)",
                          prettyNum(round(nrow(long)), big.mark = ","),
                          prettyNum(round(sum(long$mv)), big.mark = ",")), "\n",
                  sprintf("%-7s %6s %12s",
                          "net (s)",
                          prettyNum(round(nrow(short)), big.mark = ","),
                          prettyNum(round(sum(short$mv)), big.mark = ",")), "\n",
                  sprintf("%-7s %6s %12s",
                          "net   ",
                          prettyNum(round(nrow(object@actual)), big.mark = ","),
                          prettyNum(round(sum(object@actual$mv)), big.mark = ",")), "\n",
                  sprintf("%-7s %6s %12s",
                          "tot   ",
                          prettyNum(round(nrow(object@actual)), big.mark = ","),
                          prettyNum(round(sum(abs(object@actual$mv))), big.mark = ",")), "\n",
                  sep = "")
              
              
              ## For each side, display the top 10 trades by rank
              ## quality.  For buys and covers higher ranks are
              ## better; for sells and shorts, lower ranks are more
              ## attractive.
              
              cat("\nBest/worst by side/ranks quality\n\n")
              actual.sides <- split(object@actual, object@actual$side)
              for(side in all.sides[all.sides %in% names(actual.sides)]){
                
                if(side %in% c("C","B"))
                  decreasing <- TRUE
                else
                  decreasing <- FALSE

                single.side <- actual.sides[[side]]
                single.side <- single.side[order(single.side$rank.t,
                                                 decreasing = decreasing),]
                
                cat("Side: ", side, "\n", sep = "")
                num <- 5
                if(nrow(single.side) < num * 2){
                  show(single.side)
                }
                else{
                  single.side <- rbind(head(single.side, n = num),
                                       tail(single.side, n = num))
                  show(single.side)
                }
                cat("\n")
              }

              cat("\nLargest orders:\n\n")
              show(head(object@actual[order(abs(object@actual$mv), decreasing = TRUE),], n = 10))

            }
            else{
              cat("\n\nNo actual trades were generated in this tradelist.\n\n")
            }

            ## Display to the top 10 swaps by rank gain.  Row names
            ## for the swaps data frame contain symbol information
            ## of the form 'enter chunk tag,exit chunk tag".

            if (nrow(object@swaps.actual) > 0){
              swap.summary.cols <-
                c("side.enter","tca.rank.enter",
                  "side.exit","tca.rank.exit","rank.gain")
              
              swaps.actual <-
                object@swaps.actual[order(object@swaps.actual$rank.gain,
                                          decreasing = TRUE),]
              
              cat("\nBest rank gain swaps:\n\n")
              show(head(unique(subset(swaps.actual, swaps.actual$id.enter != .dummy.id() &
                                      swaps.actual$id.exit != .dummy.id(),
                                      select = swap.summary.cols)),
                        n = 5))
              
              ## Sort the other way, and look at the worst swaps.
              
              swaps.actual <-
                object@swaps.actual[order(object@swaps.actual$rank.gain,
                                          decreasing = FALSE),]

              cat("\n")
              cat("Worst rank gain swaps:\n\n")
              show(head(unique(subset(swaps.actual, swaps.actual$id.enter != .dummy.id() &
                                      swaps.actual$id.exit != .dummy.id(),
                                      select = swap.summary.cols)),
                        n = 5))
            }


            if(!is.null(object@sorts) && nrow(object@ranks)){
              cat("\n")
              cat("Sort thresholds:\n\n")

              all.sorts  <- do.call(rbind, lapply(object@rank.sorts, function(x) { x[c("id","rank")] }))
              top.ranks  <- aggregate(all.sorts[c("rank")], by = list(id = all.sorts$id), min)
              
              for(s in names(object@sorts)){

                sort.weight   <- object@sorts[[s]]
                rank.sort     <- object@rank.sorts[[s]]
                rank.sort     <- rank.sort[order(rank.sort$rank),]
                rank.sort.tmp <- merge(rank.sort, top.ranks, by = "id", all.x = TRUE)
                rank.sort.tmp <- rank.sort.tmp[order(rank.sort.tmp$rank.x),]

                ## Traded indices are the indices of orders within the
                ## sort where 1) each order (or a portion of it)
                ## appears in the 'actual' slot and 2) the order
                ## appears in no other sort with a better rank.
                
                traded.indices <- which(rank.sort.tmp$rank.x == rank.sort.tmp$rank.y &
                                        rank.sort$id %in% object@actual$id)


                rank.cutoff <- ifelse(length(traded.indices > 0), max(traded.indices), NA)
                
                cat(paste(s, ": ", rank.cutoff, sep = ""), "\n")
                
                ## If we didn't trade anything from this sort, just
                ## show the top and bottom two.

                if(length(traded.indices) == 0){
                  if(nrow(rank.sort) <= 4){
                    show(rank.sort)
                  }
                  else{
                    show(rbind(head(rank.sort, n = 2),
                               tail(rank.sort, n = 2)))
                  }
                  cat("\n")
                  next()
                }

                ## Identify the cutoff order with an asterisk.

                rank.sort$id[rank.cutoff] <- paste("(*)", rank.sort$id[rank.cutoff])

                ## If there are fewer than 4 entries in the sort, we
                ## don't have to go through the trouble of getting the
                ## best, worst, etc. and worry about indexing out of
                ## bounds.
                
                if(nrow(rank.sort) < 4){
                  show(rank.sort)
                  cat("\n")
                  next()
                }

                rank.sort.show   <- NULL
                
                ## Show the best rank if not shown already.

                if(rank.cutoff - 1 > 1){
                  rank.sort.show <- rbind(rank.sort.show, rank.sort[1,])
                }

                ## Show the cutoff point, the one better and two worse orders.
                
                rank.sort.show   <- rbind(rank.sort.show,
                                          rank.sort[(rank.cutoff - 1):(min(rank.cutoff + 2, nrow(rank.sort))),])

                ## Show the worst rank if not shown already and is not
                ## the same as the best rank.

                if(rank.cutoff + 2 < nrow(rank.sort)){
                  rank.sort.show <- rbind(rank.sort.show, rank.sort[nrow(rank.sort),])
                }
                
                show(rank.sort.show)
                cat("\n")
              }              
            }
          }
          )

## Create a data frame of dummy chunks for a given side and number.

setMethod("dummyChunks",
          "tradelist",
          function(object, side, num, quality){
            
            stopifnot(quality %in% c("good","bad"),
                      num >= 0)
            
            dummy.chunk <- data.frame(id = 0)
            dummy.chunk[chunksCols(object)] <- NA

            dummy.chunk$id   <- .dummy.id()
            dummy.chunk$side <- side
            dummy.chunk[c("orig","target","shares",
                          "mv","rank.t","tca.rank",
                          "chunk.shares","chunk.mv",
                          "chunk")] <- 0

            ## Any large number (except +/-Inf) could be used for
            ## the rank.t of dummy trades
            
            dummy.chunk$rank.t <- switch(side,
                                         B = 10000,
                                         S = -10000,
                                         C = 10000,
                                         X = -10000,
                                         stop(paste("Invalid side: ",
                                                    side, sep = "")))

            ## The reverse dummy rank should be used if bad chunks
            ## are desired.
            
            if(quality == "bad")
              dummy.chunk$rank.t <- -1 * dummy.chunk$rank.t
            else if(quality != "good")
              stop("Invalid value for quality")

            dummy.chunk$tca.rank <- dummy.chunk$rank.t

            ## Record the dummy chunk quality in the dummy.quality
            ## column.

            dummy.chunk$dummy.quality <- quality

            ## Create the correct number of dummy chunks by repeated
            ## indexing.
            
            dummy.chunk[rep(1, num),]
          }
          )

## Trade-cost-adjusted alpha

## This area of the code is messy and will change in upcoming
## versions.  Instead of supplying a list of function names to the
## 'tca' slot, the user will be able to supply a list of function
## objects that can be used for rank trade cost adjustment and won't
## be limited to the arbitrary adjustment functions included here.

## Our current naive trade code adjustment function adds an
## adjustment factor to the supplied alpha based on the 'average
## volume percentile' of the supplied chunk where

## 'average volume percentile' ==
##    mean(start of chunk as a percentage of volume,
##         end of chunk as a percentage of volume)
##

## For buys/covers, a negative number is added to alpha to make
## trades futher into a day's volume less attractive; for
## sells/short, a positive number is added to achieve the same
## affect.

## The magnitude of the additive adjustment is determined by the
## depth of the chunk within a day's volume.  Currently, we use the
## following scheme:

## Pct range    |   adjustment
## ---------------------------
## [0,    0.15) |         0.00
## [0.15, 0.30) |     +/- 0.75
## [0.30, 0.50) |     +/- 2.00
## [0.50, +Inf) |     +/- 6.00

.tca.volume <- function(alpha, shares, cumsum.shares, side, volume, chunk){

  vol.avg.pct <- mean(c((cumsum.shares - shares) / volume,
                        cumsum.shares            / volume))

  adjustment <- 0
  if((0 <= vol.avg.pct && vol.avg.pct < 0.15) || chunk == 1){
    adjustment <- 0
  }
  else if(0.15 <= vol.avg.pct && vol.avg.pct < 0.30){
    adjustment <- 0.75
  }
  else if(0.30 <= vol.avg.pct && vol.avg.pct < 0.50){
    adjustment <- 2
  }
  else if(0.50 <= vol.avg.pct){
    adjustment <- 6
  }
  else{
    stop("Encountered negative vol.avg.pct")
  }

  if(side %in% c("C","B")){
    adjustment <- -1 * adjustment
  }else if(side %in% c("X","S")){
    adjustment <- adjustment
  }
  else{
    stop(paste("Invalid side: ", side, sep = ""))
  }

  alpha + adjustment
}

## A conservative version of .tca.volume

.tca.volume.conservative <- function(alpha, shares, cumsum.shares, side, volume, chunk){

  vol.avg.pct <- mean(c((cumsum.shares - shares) / volume,
                        cumsum.shares            / volume))

  adjustment <- 0
  if((0 <= vol.avg.pct && vol.avg.pct < 0.15) || chunk == 1){
    adjustment <- 0
  }
  else if(0.15 <= vol.avg.pct && vol.avg.pct < 0.30){
    adjustment <- 3
  }
  else if(0.30 <= vol.avg.pct){
    adjustment <- 6
  }
  else{
    stop("Encountered negative vol.avg.pct")
  }

  if(side %in% c("C","B")){
    adjustment <- -1 * adjustment
  }else if(side %in% c("X","S")){
    adjustment <- adjustment
  }
  else{
    stop(paste("Invalid side: ", side, sep = ""))
  }

  alpha + adjustment
}

## Adjustment function that makes only the first 15% of the day's
## volume doable.

.tca.volume.15.pct.only <- function(alpha, shares, cumsum.shares, side, volume, chunk){

  vol.avg.pct <- mean(c((cumsum.shares - shares) / volume,
                        cumsum.shares            / volume))

  adjustment <- 0
  if((0 <= vol.avg.pct && vol.avg.pct < 0.15) || chunk == 1){
    adjustment <- 0
  }
  else{
    adjustment <- 50000
  }

  if(side %in% c("C","B")){
    adjustment <- -1 * adjustment
  }else if(side %in% c("X","S")){
    adjustment <- adjustment
  }
  else{
    stop(paste("Invalid side: ", side, sep = ""))
  }

  alpha + adjustment
}

## Another for only the first 10%.  Refactoring here would be nice.

.tca.volume.10.pct.only <- function(alpha, shares, cumsum.shares, side, volume, chunk){

  vol.avg.pct <- mean(c((cumsum.shares - shares) / volume,
                        cumsum.shares            / volume))

  adjustment <- 0
  if((0 <= vol.avg.pct && vol.avg.pct < 0.10) || chunk == 1){
    adjustment <- 0
  }
  else{

    ## The below will cause the tca.rank of all chunks above 10pct
    ## dvol to be excluded from tradelist analysis.
    
    adjustment <- NA
  }

  if(side %in% c("C","B")){
    adjustment <- -1 * adjustment
  }else if(side %in% c("X","S")){
    adjustment <- adjustment
  }
  else{
    stop(paste("Invalid side: ", side, sep = ""))
  }

  alpha + adjustment
}

## Trade-cost adjust the first chunk of new positions (positions with
## original shares == 0) in order to force their completion.

.tca.entry <- function(alpha, side, orig, chunk){
  adjustment <- 0
  if(orig == 0 && chunk == 1){
    adjustment <- 3
  }
  if(side %in% c("X","S")){
    adjustment <- -1 * adjustment
  }

  alpha + adjustment
}

.tca.too.short <- function(alpha, side, mv, mv.long.orig){
  adjustment <- 0
  if(side == "C" && (mv / mv.long.orig) > 0.0005){
    adjustment <- 3
  }
  alpha + adjustment
}

.dummy.id <- function(){
  "123456789"
}

setMethod("securityInfo",
          signature(object = "tradelist", id = "character"),
          function(object, id){
            x <- object@data

            ## Try to match on the id column first, then symbol if
            ## present in the object's data slot.

            pat <- paste("^", id, "$", sep = "")
            if(any(grep(pat, x$id))){
              y <- x[x$id == id,]
            }
            else if("symbol" %in% names(x) &&
                    any(grep(pat, x$symbol))){
              y <- x[x$symbol == id,]
            }
            else{
              stop("id not found")
            }

            info <- list()
            swaps.cols <- c("id.enter","side.enter","chunk.mv.enter","tca.rank.enter",
                            "id.exit","side.exit","chunk.mv.exit","tca.rank.exit",
                            "rank.gain")
            
            info$restricted <- object@restricted[object@restricted$id == y$id,]
            info$candidates <- object@candidates[object@candidates$id == y$id,]
            info$ranks      <- object@ranks[object@ranks$id == y$id,]
            info$chunks     <- object@chunks[object@chunks$id == y$id,]
            info$swaps      <- object@swaps[object@swaps$id.enter == y$id |
                                            object@swaps$id.exit  == y$id,]

            if(nrow(info$swaps) > 0){
              info$swaps      <- info$swaps[swaps.cols]
            }
            
            info$swaps.actual <-
              object@swaps.actual[object@swaps.actual$id.enter == y$id |
                                  object@swaps.actual$id.exit  == y$id,]

            if(nrow(info$swaps.actual) > 0){
              info$swaps.actual <- info$swaps.actual[swaps.cols]
            }
            
            info$chunks.actual <-
              object@chunks.actual[object@chunks.actual$id == y$id,]
            info$actual     <- object@actual[object@actual$id == y$id,]

            ## Display applicable sort information
            
            info$rank.sorts <- lapply(object@rank.sorts, function(x) { x[x$id == y$id,] })
            info$rank.sorts <- info$rank.sorts[sapply(info$rank.sorts, nrow) > 0]
            
            show(info)
            invisible(info)
            
          }
          )
