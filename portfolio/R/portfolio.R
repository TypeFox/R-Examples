################################################################################
##
## $Id: portfolio.R 1311 2008-10-31 17:38:03Z enos $
##
## A more complex, full-featured portfolio object that includes
## shares, a notion of equity, and is better suited for use in
## implementation.
##
################################################################################

setMethod("initialize",
          signature(.Object = "portfolio"),
          function(.Object, ...){
            .Object <- callNextMethod()

            if(nrow(.Object@shares) == 0){
              row.names(.Object@shares) <- integer(0)
            }
            .Object
          }
          )

## Calculate weights from shares.
##
## Calculate a weight for each position and store in the 'weights'
## slot data frame.  Weight of a position is calculated with respect
## to the market value of all positions on the same side as the
## position.

setMethod("calcWeights",
          "portfolio",
          function (object){

            ## There are some checks on weight.style and equity in the
            ## validity method for class 'portfolio', but we repeat
            ## some here to be safe.  This method is usually called on
            ## an object that isn't valid (after, say, specifying
            ## shares and wanting to calculate weights to make the
            ## object valid).  Is there a better way to handle this?

            stopifnot(is.character(object@weight.style) &&
                      length(object@weight.style) == 1 &&
                      object@weight.style %in% c("sides.separate", "long.tmv",
                                                 "short.tmv", "reference.equity"))
            
            ## Nothing to do here if there are no positions.
            
            if(nrow(object@shares) == 0)
              return(object)

            x <- merge(object@shares, object@data, by = "id", all.x = TRUE)
            
            ## Can't calculate mv without price.var.

            if(! object@price.var %in% names(x)){
              stop("Cannot calculate weights without a valid price.var!")
            }

            x$mv <- x$shares * x[[object@price.var]]

            mv.long  <- mvLong(object)
            mv.short <- mvShort(object)

            if(isTRUE(all.equal(object@weight.style, "sides.separate"))){
              weight.denom <- ifelse(x$shares > 0, mv.long, abs(mv.short))
            }
            else if(isTRUE(all.equal(object@weight.style, "long.tmv"))){
              weight.denom <- mv.long
            }
            else if(isTRUE(all.equal(object@weight.style, "short.tmv"))){
              weight.denom <- abs(mv.short)
            }
            else if(isTRUE(all.equal(object@weight.style, "reference.equity"))){

              if(is.null(object@equity) ||
                 length(object@equity) != 1){
                stop(paste("Valid equity slot contents required",
                           "for reference.equity weight.style"))
              }
              
              weight.denom <- object@equity
            }
            else{
              stop("Invalid weight.style")
            }
            
            x$weight <- x$mv / weight.denom
            
            object@weights <- x[c("id","weight")]

            ## As a final check, ensure that the weights and shares
            ## slots contain exactly the same set of securities.

            xx <- merge(object@weights, object@shares, by = "id")
            stopifnot(nrow(xx) == nrow(object@shares))
            
            invisible(object)
          }
          )

## Calculate shares from weights.
##
## Shares are calculated from weights by multiplying portfolio equity
## by weight for each security.  Additional adjustments to share
## values can be made if the necessary data is present, such as
## adjusting to round lot sizes.

setMethod("calcShares",
          "portfolio",
          function(object){
            
            ## Nothing to do here if there are no positions.
            
            if(nrow(object@weights) == 0)
              return(object)

            if(is.na(object@equity) || object@equity <= 0){
              stop("Must have non-zero positive equity to calculate shares!")
            }

            x <- merge(object@weights, object@data, by = "id", all.x = TRUE)
            
            ## Can't calculate shares without price.var.

            if(! object@price.var %in% names(x)){
              stop("Cannot calculate shares without a valid price.var!")
            }

            x$mv     <- x$weight * object@equity
            x$shares <- x$mv / x[[object@price.var]]

            ## If round.lot is not present, set all round lots to 1.
            ## We need to add a round lot adjustment on/off switch
            ## somewhere.

            if(nrow(x) > 0 && !"round.lot" %in% names(x)){
              x$round.lot <- 1
            }
            x$shares <- .nearest.multiple(x$shares, x$round.lot)
            
            object@shares <- x[c("id","shares")]

            ## As a final check, ensure that the weights and shares
            ## slots contain exactly the same set of securities.

            xx <- merge(object@weights, object@shares, by = "id")
            stopifnot(nrow(xx) == nrow(object@weights))
            
            invisible(object)
          }
          )

setMethod("create",
          "portfolio",
          function(object, ...){

            object <- callNextMethod()
            object <- calcShares(object)

            object
          }
          )

## Long/Short market value calculations based on shares and price.
## Note that a valid object is not required here, to aid such
## situations as calculating weights for the first time from a set of
## shares.

setMethod("mvLong",
          "portfolio",
          function(object){
            
            if(nrow(object@shares) == 0)
              return(0)

            x <- merge(object@data, object@shares, by = "id")
            
            ## Can't calculate mv without price.var.

            if(! object@price.var %in% names(x)){
              stop("Cannot calculate market value without a valid price.var!")
            }

            x$mv <- x$shares * x[[object@price.var]]
            sum(x$mv[x$mv > 0], na.rm = TRUE)
          }
          )

setMethod("mvShort",
          "portfolio",
          function(object){
            
            if(nrow(object@shares) == 0)
              return(0)

            x <- merge(object@data, object@shares, by = "id")
            
            ## Can't calculate mv without price.var.

            if(! object@price.var %in% names(x)){
              stop("Cannot calculate market value without a valid price.var!")
            }

            x$mv <- x$shares * x[[object@price.var]]
            sum(x$mv[x$mv < 0], na.rm = TRUE)
          }
          )

## Currently I don't provide these methods for the portfolioBasic
## class, although I probably should.  I also should think about a
## better way of dealing with sides here.  For instance, if you have a
## long-onlyportfolio it's annoying to have to call sizeLong and
## mvLong the whole time.

setMethod("sizeLong",
          "portfolio",
          function(object){
            return(length(object@shares$shares[object@shares$shares > 0]))
          }
          )

setMethod("sizeShort",
          "portfolio",
          function(object){
            return(length(object@shares$shares[object@shares$shares < 0]))
          }
          )

## Try to look up information for this security in this portfolio
## using "x" as a lookup identifier.

setMethod("securityInfo",
          signature(object = "portfolio", id = "character"),
          function(object, id){

            x <- merge(object@data, object@weights, by = "id")
            x <- merge(x, object@shares)

            ## Try to match on the id column first, then symbol if
            ## present in the object's data slot.
            
            pat <- paste("^", id, "$", sep = "")           
            if(any(grep(pat, x$id))){
              y <- x[x$id == id,]
            }
            else if(length(object@symbol.var) > 0 &&
                    object@symbol.var %in% names(x) &&
                    any(grep(pat, x[[object@symbol.var]]))){
              y <- x[x[[object@symbol.var]] == id,]
            }
            else{
              stop("id not found")
            }

            if(length(object@price.var) == 1 &&
               object@price.var %in% names(x)){
              y$mv <- y[[object@price.var]] * y$shares
            }

            columns <- c("id", "weight", "shares",
                         "mv", object@in.var)

            if(length(object@symbol.var) > 0 &&
               object@symbol.var %in% names(x)){
              row.names(y) <- y[[object@symbol.var]]
            }

            
            show(y[columns[columns %in% names(y)]])
            
            invisible(y)
          }
          )

setMethod("all.equal",
          signature(target = "portfolio", current = "portfolio"),
          function(target, current){
            
            r <- callNextMethod()

            if(!isTRUE(r)){
              return(r)
            }

            s1 <- target@shares
            s2 <- current@shares

            s1 <- s1[order(s1$id),]
            s2 <- s2[order(s2$id),]

            if(!isTRUE(all.equal(s1$id, s2$id))){
              return(paste("Identifiers:", all.equal(s1$id, s2$id)))
            }
            else if(!isTRUE(all.equal(s1$shares, s2$shares))){
              return(paste("Shares:", all.equal(s1$shares, s2$shares)))
            }
            else TRUE
          }
          )

setMethod("+",
          signature(e1 = "portfolio", e2 = "portfolio"),
          function(e1, e2){

            validObject(e1)
            validObject(e2)
            
            r.basic <- callNextMethod()
            r <- new("portfolio", type = "unknown", size = "unknown")
            r@data <- r.basic@data
            r@weights <- r.basic@weights

            r@id.var <- r.basic@id.var
            r@symbol.var <- r.basic@symbol.var
            r@ret.var <- r.basic@ret.var

            if(nrow(e1@shares) == 0 && nrow(e2@shares) == 0){
              return(r)
            }
            else if(nrow(e1@shares) == 0){
              r@shares <- e2@shares
              return(r)
            }
            else if(nrow(e2@shares) == 0){
              r@shares <- e1@shares
              return(r)
            }
            else{
              
              w <- merge(subset(e1@shares, !is.na(e1@shares$shares)),
                         subset(e2@shares, !is.na(e2@shares$shares)),
                         suffixes = c(".e1", ".e2"), by = "id", all = TRUE)
              

              w$shares.e1[is.na(w$shares.e1)] <- 0
              w$shares.e2[is.na(w$shares.e2)] <- 0
              
              w$shares <- w$shares.e1 + w$shares.e2
              r@shares <- w[c("id","shares")]

              ## Remove entries that now have zero shares.

              r@shares <- r@shares[is.na(r@shares$shares) | r@shares$shares != 0,]
            }
            return(r)
          }
          )

setMethod("getYahooData",
          signature(object = "portfolio", symbol.var = "character"),
          function(object, symbol.var = object@symbol.var){

            data <- object@data

            ## For data.frame merges to work correctly, all symbol
            ## names must be unique

            if(length(unique(data[[symbol.var]])) != nrow(data)){
              warning("Not all symbols names are unique.")
              return(object)
            }
            
            ## The symbols may not contain "@", "/", "+", "?", or "-"
            ## characters, otherwise Yahoo! will not interpret the request
            ## string correctly

            invalid.names <- grep("[@/+&?-]", data[[symbol.var]])

            if(any(invalid.names)){

              ## subsets out the invalid names and warns the user

              data <- data[-invalid.names,]

              warning(
                      cat("Symbols contain invalid characters. Not retrieving data for these stocks: \n",
                          paste(data[grep("[@/+&?-]", data[[symbol.var]]), symbol.var], sep = ", ")
                          )
                      )

            }

            ##  builds query string as documented in Perl module Finance:YahooQuote

            QURL.base <- "http://quote.yahoo.com/d?f="
            QURL.format <- "srb4j1p6p5"
            QURL.end <- "&s="

            QURL.format <- paste(QURL.base, QURL.format, QURL.end, sep = "")

            ## retrieves data according to the stock's symbol

            my.symbols <- object@data[[symbol.var]]
            my.symbols <- as.character(my.symbols)

            ## Creates a temporary file where we store the data
            ## downloaded from Yahoo!

            my.csv <- tempfile()
            on.exit(unlink(my.csv))
            
            while(length(my.symbols) > 0){

              ## Yahoo onlys allows the request of 199 stocks at a time

              if(length(my.symbols) > 199){

                sub.symbols <- my.symbols[1:199]
                my.symbols <- my.symbols[200:length(my.symbols)]

              }else{

                sub.symbols <- my.symbols
                my.symbols <- my.symbols[0]

              }

              ## builds the query string

              q.string <- paste(sub.symbols, collapse = "+")
              QURL <- paste(QURL.format, q.string, sep = "")

              ## uses "append" mode in case we have more than 200 symbols

              status <- download.file(QURL, destfile = my.csv,
                                      method = "auto", mode = "a", quiet = TRUE)

              if(status != 0){

                stop(paste("download error, status", status))
              }  
            }
            
            ## interprets the data

            col.names <- c(symbol.var, "p.e.ratio", "book.value",
                           "market.cap.bil", "price/book","price/sales")

            colClasses <- c("character", "numeric", "numeric", "character",
                            "numeric", "numeric")

            x <- read.table(my.csv, sep = ",", col.names = col.names, as.is = TRUE,
                            blank.lines.skip = TRUE, comment.char = "", nrows = nrow(data),
                            na.strings = "N/A")
            
            ## Removes the "B" from the market.cap.bil column and converts the
            ## values to numeric

            x[["market.cap.bil"]] <- as.character(x[["market.cap.bil"]])

            x[["market.cap.bil"]] <- gsub('B$', '', x[["market.cap.bil"]])
            
            x[["market.cap.bil"]] <- as.numeric(x[["market.cap.bil"]])

            object@data <- merge(object@data, x, by = symbol.var)
            
            object
          }

          )

## "Expose" a portfolio with a data.frame of trades, that is, apply a
## set of trades to the set of positions (shares) of an existing
## portfolio and return the resulting portfolio.

## It's pretty clear that we shouldn't automatically recalculate
## weights in this method.  However, I believe we should return an
## object that is valid.  Invalidity can arise when trades close or
## open a position.  After that occurs, the weights data frame will
## have a different set of securities than the shares data frame.

## For now, we add a row to the weights data frame for newly opened
## positions with NA weight.  Positions that are closed by expose()
## have their corresponding weights data frame entries removed.  All
## other values in the weights data frame stay the same.

setMethod("expose",
          signature(object = "portfolio", trades = "trades"),
          function(object, trades){
            
            ## Verify we're working with valid objects:

            validObject(object)
            validObject(trades)
            
            trades <- trades@trades
            
            if(nrow(trades) == 0){
              return(object)
            }
                        
            ## The trades data frame can contain side changes, which
            ## complicates matters here a bit.  First sort by id,
            ## putting side exits first.  Then collect the entries so
            ## that we can apply later.
            
            trades    <- trades[order(trades$id, match(trades$side, c("S","C","B","X"))),]
            trades.sc <- trades[duplicated(trades$id),]
            trades    <- trades[!duplicated(trades$id),]
            
            ## Here's a dot function that doesn't have to deal with
            ## multiple orders per stock, and therefore can be much
            ## more conservative about its input.  We'll need to apply
            ## this twice.
            
            .expose <- function(object, trades){
              stopifnot(all(!duplicated(trades$id)))
              
              shares <- object@shares
              
              if(nrow(shares) == 0){
                shares <- data.frame(id = I(as.character(NA)), shares = as.numeric(NA))
              }
              
              x <- merge(shares, trades, by = "id", all = TRUE,
                         suffixes = c(".orig", ".exp"))
            
              ## Just in case there are any rows with NA id (like the
              ## dummy row we had to create in the case of an empty
              ## target), remove rows with NA id.
              
              x <- subset(x, !is.na(id))

              x$shares.orig <- ifelse(is.na(x$shares.orig), 0, x$shares.orig)
              x$shares.exp  <- ifelse(is.na(x$shares.exp), 0, x$shares.exp)

              ## Now check to see if anything nasty is happening.

              if(any(x$side %in% c("B","S") & x$shares.orig < 0) ||
                 any(x$side %in% c("C","X") & x$shares.orig > 0) ||
                 any(x$side %in% c("C","S") & abs(x$shares.exp) > abs(x$shares.orig))
                 ){
                
                ## This error message can be much more informative.
                
                stop("Illegal trades found")
                
              }

              ## At this point, we're ready to perform a simple sum.

              x$shares <- x$shares.orig + ifelse(x$side %in% c("S","X"), -1, 1) * x$shares.exp

              ## There should be no NA's at this point.

              stopifnot(all(!is.na(x$shares)))

              ## Take away closed positions.

              x <- subset(x, shares != 0)

              object@shares <- x[c("id","shares")]

              ## Now clean up weights data frame.  First, remove
              ## weights for closed positions.

              object@weights <- subset(object@weights, id %in% object@shares$id)

              ## Then add NA weights for opened positions.
              
              new.weights <- subset(object@shares, !id %in% object@weights$id)
              if(nrow(new.weights) > 0){
                names(new.weights) <- c("id","weight")
                new.weights$weight <- NA
                object@weights <- rbind(object@weights, new.weights)
              }
              
              object
            }

            object <- .expose(object, trades)
            if(nrow(trades.sc) > 0){
              object <- .expose(object, trades.sc)
            }

            ## Require the caller to calc weights.
            
            ## object <- calcWeights(object)
            
            invisible(object)
          }
          )

setMethod("updatePrices",
          signature(object = "portfolio",
                    id     = "character",
                    price  = "numeric"),
          function(object, id, price){

            object@data[[object@price.var]] <-
              price[match(object@data$id, id)]

            invisible(object)
          }
          )


setMethod("performance",
          signature(object = "portfolio"),
          function(object,
                   market.data = NULL){

            ## Verify that we're working with a valid object.

            validObject(object)

            if(is.null(market.data)){
              
              ## If we're not supplied market data, use
              ## portfolioBasic's performance calculation.
              
              return(invisible(callNextMethod()))
              
            }

            perf <- new("performance")
            
            ## Performance is calculated over a period of time, using
            ## the information contained in the market.data data
            ## frame.  It is required that this data frame have the
            ## following columns: start.price, ret, end.price.

            ## At the start of the period we have a price
            ## (start.price) for each security that we use to compute
            ## each security's market value and weight.  During the
            ## period each security has a total return (ret), using
            ## which we can calculate the total profit in each stock.
            ## Finally, at the end of the period we have a price that
            ## we use to reflect market values/weights at the end of
            ## the performance period.

            ## The distiction between these two prices is useful; it's
            ## also useful to centralize the handling of prices as
            ## opposed to directly setting the price.var column of the
            ## portfolio object's data slot.

            ## The distinction is also in the right direction of
            ## requiring the full set of information required to
            ## reflect the true state of the portfolio at the end of
            ## this period.  Since we require total return, cash
            ## dividends are handled correctly.  Profit in a stock
            ## that has split will be handled correctly only if
            ## adjusted prices are supplied (which is recommended at
            ## this point).  Stock dividends, however, pose an issue
            ## because the end-of-period shares will have changed but
            ## is not handled at all in this method.  We will provide
            ## places to input such data in the future.

            ## Note that we deal in a single reference currency
            ## (unspecified) in the current version.

            ## Also note that we ignore all notion of price already
            ## contained in the portfolio object in this method.

            ## At this stage should we force client to pass in column
            ## 'id', or use id.var from the object?  Using the former
            ## for now because we may do away with id.var.
            
            stopifnot(is.data.frame(market.data),
                      all(c("id","start.price","end.price","ret")
                          %in% names(market.data)))

            ## If this portfolio is empty, return the default
            ## performance object with itself as t-plus-one.
            
            if(nrow(object@shares) == 0){
              perf@t.plus.one <- object
              return(perf)
            }

            ## First, update the copy of this portfolio (object) with
            ## the price supplied in market.data.  Perhaps I should
            ## also pass in a price.var here, but for now use the
            ## price.var specified in the object.

            object <- updatePrices(object,
                                   market.data$id,
                                   market.data$start.price)

            ## Calculate weights based on these prices.

            object <- calcWeights(object)

            ## Now collect all this information in 'x' so that we can
            ## calculate contributions and profits.
            
            x <- merge(object@weights, object@shares,
                       by = "id", all = TRUE)
            x <- merge(x, market.data, by = "id", all.x = TRUE)

            
            x$contrib     <- x$weight * x$ret

            x$mv          <- x$shares * x$start.price
            x$profit      <- x$mv * x$ret

            long.denom    <- mvLong(object)
            short.denom   <- abs(mvShort(object))
            x$profit.contrib  <- x$profit / ifelse(sign(x$mv) > 0,
                                                   long.denom,
                                                   short.denom)

            ## By construction, x$contrib and x$profit.contrib should
            ## be the name number, so this is a good sanity check that
            ## other parts of the portfolio package are functioning
            ## properly.  Note that we must be clear that the
            ## resulting return in this object assumes we're
            ## calculating weight relative to one's side.

            stopifnot(all.equal(x$contrib, x$profit.contrib))

            ## Collect results in an object of class
            ## 'performance'.

            perf@ret             <- sum(x$contrib, na.rm = TRUE)
            perf@profit          <- sum(x$profit,  na.rm = TRUE)
            perf@ret.detail      <- x[c("id", "weight", "ret", "contrib","profit")]

            perf@missing.price   <- sum(is.na(x$start.price))
            perf@missing.return  <- sum(is.na(x$ret))

            if(length(object@symbol.var) > 0 &&
               object@symbol.var %in% names(object@data)){

              symbols <- object@data[[object@symbol.var]]

              if(any(is.na(symbols[match(perf@ret.detail$id, object@data$id)]))){
                warning(paste("Attempt to use symbol.var column",
                              "to set row.names failed: NA's exist"))
              }
              else{
                row.names(perf@ret.detail) <-
                  symbols[match(perf@ret.detail$id, object@data$id)]
              }
            }
            
            ## Now use end prices to reflect the state of the
            ## portfolio at the end of the performance period.

            object <- updatePrices(object, market.data$id, market.data$end.price)
            object <- calcWeights(object)
            
            perf@t.plus.one <- object
            
            invisible(perf)
          }
          )

setMethod("contribution",
          signature(object = "portfolio", contrib.var = "character"),
          function(object, contrib.var, buckets = 5, market.data){

            validObject(object)
            
            x <- object@data

            stopifnot(all(contrib.var %in% names(x)))

            ## Call the performance method and merge results into our
            ## data.frame. 

            ## I should probably add a check here to see if there is
            ## anything in the object 'perf'.

            perf <- performance(object, market.data = market.data)
            x <- merge(x, perf@ret.detail, by = "id")
            x$weight <- abs(x$weight)
            
            result <- list()
            
            for(att in contrib.var){

              att.cut <- att
              
              if(is.numeric(x[[att]])){
                att.cut <- paste(att, "levels", sep = ".")

                ## Note that numeric category quantiles are computed
                ## over the universe, or securities in the 'data'
                ## slot, not just the values for contrib.var among
                ## this portfolio's positions.

                quantiles <- quantile(object@data[[att]],
                                      seq(0, 1, 1/buckets), na.rm = TRUE)
                
                x[[att.cut]] <- cut(x[[att]], quantiles, na.rm = TRUE)
              }
              a <- aggregate(list(x[c("weight","contrib")]),
                             list(variable = x[[att.cut]]), sum, na.rm = TRUE)
              
              ## Make sure that all levels in the universe-based
              ## intervals are present in this aggregate, and that
              ## they appear in the correct order.

              all.levels <- levels(x[[att.cut]])

              if(any(! all.levels %in% a$variable)){
                a <- rbind(a, data.frame(variable = all.levels[! all.levels %in% a$variable],
                                         weight = 0, contrib = 0))
                a <- a[match(all.levels, a$variable),]
              }
              
              ## To compute roic, return on invested capital, divide
              ## contribution by weight scaled to [0,1].

              ## ** Should I leave roic as NaN (instead of 0) for
              ## those categories with 0 weight?

              a$weight   <- a$weight / sum(a$weight)
              a$roic     <- ifelse(a$weight == 0, 0, a$contrib / a$weight)

              if(is.numeric(x[[att]])){

                a$rank <- 1:nrow(a)

                ## Don't use high/low if only one category.
                
                if(nrow(a) != 1){
                  a$rank[1] <- "1 - low"
                  a$rank[nrow(a)] <- paste(nrow(a), "high", sep = " - ")
                }

                ## Put rank first.
                
                a <- a[c(which("rank" == names(a)),
                         which("rank" != names(a)))]

              }

              result[[att]] <- a
            }

            contrib.obj <- new("contribution", data = result)
            contrib.obj
          }
          )

setMethod("portfolioDiff",
          signature(object = "portfolio", x = "portfolio"),
          function(object, x){

            ## For this method on portfolioBasic objects we require
            ## that both objects have data object with the same
            ## columns.  Here, we take the opposite approach (and
            ## perhaps we should do the same for portfolioBasic).  The
            ## resulting diff portfolio will have only columns found
            ## in both data objects.

            keep.cols <- intersect(names(object@data), names(x@data))

            ## But, keep the handling of different rows the same (the
            ## left operand's data slot gets precedence).
            
            p.diff.data <- rbind(object@data[keep.cols],
                                 subset(x@data[keep.cols], ! id %in% object@data$id))

            p.diff <- new("portfolio",
                          name = "Portfolio diff", data = p.diff.data)

            ## Again, left operand gets precedence.

            p.diff@price.var <- object@price.var

            if(nrow(object@shares) > 0 && nrow(x@shares) == 0){
              p.diff@shares <- object@shares
            }
            else if(nrow(object@shares) == 0 && nrow(x@shares) > 0){
              p.diff@shares <- x@shares
            }
            else if(nrow(object@shares) > 0 && nrow(x@shares) > 0){

              s.diff    <- merge(object@shares, x@shares, by = "id", all = TRUE)
              s.diff.na <- NULL
              
              ## I don't think it should be possible to have NA shares.
              ## Require this of the two operands, but we might want to
              ## put in validity method later.

              stopifnot(all(!is.na(object@shares$shares)),
                        all(!is.na(x@shares$shares)))
              
              ## Make NA weights 0 for those stocks that were in one
              ## portfolio but not in the other.

              if(any(is.na(s.diff$shares.x))){
                s.diff$shares.x[is.na(s.diff$shares.x)] <- 0
              }
              if(any(is.na(s.diff$shares.y))){
                s.diff$shares.y[is.na(s.diff$shares.y)] <- 0
              }
              
              s.diff$shares <- s.diff$shares.y - s.diff$shares.x

              ## Those securities that have the same shares in both
              ## portfolios should not be included in the diff
              ## portfolio.
              
              s.diff <- subset(s.diff, shares != 0)
              
              p.diff@shares <- rbind(s.diff[c("id","shares")], s.diff.na)
            }

            p.diff <- calcWeights(p.diff)
            invisible(p.diff)
            
          }
          )


## Expand the contents of the data slot to include empty records for
## all id's in the shares slot but not in the data slot.

setMethod("expandData",
          signature(object = "portfolio"),
          function(object){

            ## Let's make no claim about the shares and weights slots
            ## being in sync, but focus on the shares slot since it's
            ## the one that cannot harbor NA's.

            ids.nok <- object@shares$id[!object@shares$id %in% object@data$id]

            if(length(ids.nok) > 0){
              data.template <- object@data[FALSE,]
              data.addl     <- data.template[1:length(ids.nok),]
              data.addl$id  <- ids.nok
              object@data   <- rbind(object@data, data.addl)
            }
            object
          }
          )


setMethod("summary",
          signature(object = "portfolio"),
          function(object){

            validObject(object)
            
            weight.col <- ifelse(nrow(object@weights) > 200, "bps", "pct")
            x <- .summary.prep.df(object, weight.col)

            x$shares <- object@shares$shares[match(x$id, object@shares$id)]
            x$value  <- x$shares * x[[object@price.var]]
            
            columns  <- c("id", weight.col)
            disp.num <- 5
            
            long  <- subset(x, weight > 0)
            short <- subset(x, weight < 0)

            cat("Portfolio: ", object@name, "\n\n",
                sprintf("       %6s %12s %15s", "count", "weight", "value"), "\n",

                ifelse("long" %in% object@sides,
                       paste(
                             sprintf("Long:  %6s %12s %15s",
                                     prettyNum(nrow(long), big.mark = ","),
                                     prettyNum(sum(long$weight), big.mark = ","),
                                     prettyNum(round(sum(long$value, na.rm = TRUE)), big.mark = ",")
                                     ), "\n"),
                       ""),
                ifelse("short" %in% object@sides,
                       paste(
                             sprintf("Short: %6s %12s %15s",
                                     prettyNum(nrow(short), big.mark = ","),
                                     prettyNum(sum(short$weight), big.mark = ","),
                                     prettyNum(round(sum(short$value, na.rm = TRUE)), big.mark = ",")
                                     ), "\n"),
                       ""),
                "\n",
                "Top/bottom positions by weight:\n",
                sep = "")

            ## The following section will be expanded to include
            ## more top/bottom summaries akin to our perl portfolio
            ## summary.  Currently, the only data frame we show is
            ## sorted by weight.

            by.weight <- x[order(-x$weight, na.last = NA),]
            if(disp.num * 2 > nrow(by.weight)){
              show(by.weight[columns])
            }
            else{
              show(rbind(head(by.weight, n = disp.num),
                         tail(by.weight, n = disp.num))[columns])
            }
            
            if(nrow(x[is.na(x$weight),]) > 0){
              cat("\nNA weight:\n")
              show(x[is.na(x$weight),][columns])
            }
          }
          )
