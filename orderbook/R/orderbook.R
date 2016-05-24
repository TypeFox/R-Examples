################################################################################
##
##
## orderbook.R: functions of the orderbook object
##
## limitob is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## limitob is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with limitob.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

setClass("orderbook", representation(current.ob   = "data.frame",
                                     current.time = "numeric",
                                     file         = "character",
                                     file.index   = "numeric",
                                     ob.data      = "hash",
                                     trade.data   = "vector",
                                     trader       = "logical"
                                     ),

         prototype(current.ob   = data.frame(),
                   current.time = 0,
                   file		= character(),
                   file.index   = 1,
                   ob.data      = hash(),
                   trade.data   = vector(),
                   trader       = FALSE
                   )
         )

## Reads the next n orders from the file. Use negative n to go
## backwards, but this really just reads everything over again.

setMethod("read.orders",
          signature(object = "orderbook"),
          function(object, n = 1000){
              if(identical(n, 0))
                  invisible(object)
              else if(n > 0){
                  invisible(.read.orders(object, n))
              } else {
                  n <- object@file.index + n - 1
                  object <- reset(object)
                  invisible(.read.orders(object, n))
              }

          }
          )

## Reads orders from the file until the time specified.

setMethod("read.time",
           signature(object = "orderbook"),
           function(object, n){

               ## If the time you are reading is greater than current
               ## time there is no reason to start from the beginning.
               ## This takes care of that.

               if(.to.ms(n) > object@current.time){

                   ## get.time.row finds the row in the file.

                   n <- .get.time.row(object@file, .to.ms(n),
                                      object@file.index)

                   n <- n - object@file.index

                   invisible(read.orders(object, n))
               } else {
                   n <- .get.time.row(object@file, .to.ms(n))
                   object <- reset(object)
                   invisible(read.orders(object, n))
               }
           }
           )




## Returns a vector with the price and size of the best bid order at top
## priority.

setMethod("best.bid",
          signature(object = "orderbook"),
          function(object, ...){
              x <- object@current.ob

              ## Takes out the bids.

              x <- x[x[["type"]] == "BID",]

              ## Sorts by time.

              x <- x[order(x[["time"]]),]

              ## Gets indices for the best bid.

              index <- x[["price"]] == max(x[["price"]])

              price <- x[["price"]][index]
              size <- x[["size"]][index]

              ## Return -1 (no best bid) if x is empty, otherwise return named
              ## vector of price and size.

              if(identical(nrow(x), 0)){
                  return(-1)
              } else {
                  return(c(price = price[1], size = size[1]))
              }
          }
          )

## Returns a vector with the price and size of the best ask order at top
## priority. See above for comments

setMethod("best.ask",
          signature(object = "orderbook"),
          function(object, ...){
              x <- object@current.ob

              x <- x[x[["type"]] == "ASK",]

              x <- x[order(x[["time"]]),]

              index <- x[["price"]] == min(x[["price"]])

              price <- x[["price"]][index]
              size <- x[["size"]][index]

              if(identical(nrow(x), 0)){
                  return(-1)
              } else {
                  return(c(price = price[1], size = size[1]))
              }

          }
          )

## Show basic information about the order book.

setMethod("show",
          signature(object = "orderbook"),
          function(object){

              if(isTRUE(object@trader))
                  cat("An object of class orderbook (trader)\n")
              else
                  cat("An object of class orderbook (default)\n")

              cat("--------------------------\n")
              cat("Current orderbook time:   ",
                  .to.time(object@current.time), "\n")
              cat("Message Index:            ",
                  .prettify(object@file.index - 1, "s"), "\n")
              cat("Bid Orders:               ",
                  .prettify(bid.orders(object), "s"), "\n")
              cat("Ask Orders:               ",
                  .prettify(ask.orders(object), "s"), "\n")
              cat("Total Orders:             ",
                  .prettify(total.orders(object), "s"), "\n")
          }
          )


## Plot, basically just calls the plot method. See orderbook.plot.R.

setMethod("plot",
          signature(x = "orderbook"),
          function(x, bounds = 0.1, n = 10, type = "n"){
              if(isTRUE(type %in% "n")){

                  tmp <- .plot.ob(x, bounds)
                  print(tmp)

              } else if(isTRUE(type %in% "s")){

                  .plot.side.ob(x, n)

              } else if(isTRUE(type %in% "o")){

                  tmp <- .plot.orders.ob(x, bounds)
                  print(tmp)

              } else if(isTRUE(type %in% "sd")){

                  tmp <- .supply.demand.plot(x, bounds)
                  print(tmp)

              } else if(isTRUE(type %in% "t")){

                  tmp <- .animate.plot(x, bounds)
                  print(tmp)

              } else {

                  print("Invalid type")

              }
          }
          )

## Displays summary information.

setMethod("summary",
          signature(object = "orderbook"),
          function(object){
              cat("\nCurrent time is",
                      .to.time(object@current.time), "\n\n")
              cat("Ask price levels:  ",
                  .prettify(ask.price.levels(object), "s"), "\n")
              cat("Bid price levels:  ",
                  .prettify(bid.price.levels(object), "s"), "\n")
              cat("Total price levels:",
                  .prettify(total.price.levels(object), "s"), "\n")
              cat("-----------------------------\n")
              cat("Ask orders:        ",
                  .prettify(ask.orders(object), "s"), "\n")
              cat("Bid orders:        ",
                  .prettify(bid.orders(object), "s"), "\n")
              cat("Total orders:      ",
                  .prettify(total.orders(object), "s"), "\n")
              cat("-----------------------------\n")
              cat("Spread:            ",
                  .prettify(spread(object)), "\n\n")

              ## Midpoint--format it if there are less than 3 decimals

              mid <- mid.point(object)
              check <- mid - floor(mid)
              if(!identical(nchar(check), 5))
                 formatC(mid, format = "f", digits = 2)

              cat("Mid point:         ", mid, "\n")
              cat("-----------------------------\n")
              cat("Inside market \n \n")

              ## Calls inside.market function which prints it

              inside.market(object)
              cat("\n")
          }
          )

## Displays the price levels and sizes. n specifies the number of rows to be
## displayed for ask and bid. Short = FALSE returns the data frame sorted
## in decreasing order by price so you can see all individual orders.

setMethod("display",
          signature(object = "orderbook"),
          function(object, n = 5, short = TRUE, ...){
              if(short){
                  x <- .combine.size(object, Inf)

                  ## Create ask and bid data frames

                  ask <- x[x[["type"]] == "ASK",]
                  bid <- x[x[["type"]] == "BID",]

                  cat("\nCurrent time is",
                      .to.time(object@current.time), "\n\n")
                  cat("\t\t Price \t Ask Size\n")
                  cat("---------------------------------------------\n")

                  ## Print out the top n ask prices/sizes

                  for(i in rev(1:min(n, nrow(ask)))){
                      cat("\t\t",
                          .prettify(ask[["price"]][i]), "\t",
                          .prettify(ask[["size"]][i], "s"), "\n")
                  }
                  cat("---------------------------------------------\n")

                  ## Print out the top n bid prices/sizes

                  for(i in rev(max(1, nrow(bid) - n + 1):nrow(bid))){

                      ## Spacing for right alignment, max size is
                      ## 10 million for this to work

                      size = .prettify(bid[["size"]][i], "s")
                      if(nchar(size) <= 7){
                          space = rep(" ", 7 - nchar(size))
                      } else {
                          space = ""
                      }

                      space = paste(space, collapse = "")

                      size = paste(space, size, sep = "")

                      ## Actually printing it out

                      cat(size, "\t",
                          .prettify(bid[["price"]][i]), "\n")
                  }
                  cat("---------------------------------------------\n")
                  cat("Bid Size \t Price\n")
                  invisible(object)
              } else {

                  ## Print the whole thing short == FALSE

                  x <- object@current.ob
                  x <- x[order(x[["price"]], decreasing = TRUE),]
                  return(x)
              }


          }
          )

## Go to next trade after current time.

setMethod("next.trade",
          signature(object = "orderbook"),
          function(object){
              ## .get.next.trade finds the next trade

              n <- .get.next.trade(object@file, object@file.index)
              n <- n - object@file.index
              invisible(read.orders(object, n))
          }
          )

## Go to the first trade to occur before the current time.

setMethod("previous.trade",
          signature(object = "orderbook"),
          function(object){

              if(isTRUE(object@trader))
                  skip <- 5
              else
                  skip <- 4

              x <- object@trade.data

              trdrow = as.integer(x[length(x) - skip])

              if(trdrow >= object@file.index)
                  trdrow = as.integer(x[length(x) - skip * 2 - 1])

              invisible(read.orders(object, trdrow - object@file.index))

          }
          )


## Returns the number of bid price levels.

setMethod("bid.price.levels",
          signature(object = "orderbook"),
          function(object, ...) {
              x <- object@current.ob

              ## Pull out the bids, then split by price.

              x <- x[x[["type"]]=="BID",]

              ## Use table to count.

              return(length(table(x[["price"]], exclude = NA)))
          }
          )

## Returns the number of ask price levels.

setMethod("ask.price.levels",
          signature(object = "orderbook"),
          function(object, ...) {
              x <- object@current.ob

              x <- x[x[["type"]]=="ASK",]

              return(length(table(x[["price"]], exclude = NA)))

          }
          )

## Returns the total number of price levels.

setMethod("total.price.levels",
          signature(object = "orderbook"),
          function(object, ...) {

              return(bid.price.levels(object) +
                     ask.price.levels(object))

          }
          )

## Returns the number of bid orders.

setMethod("bid.orders",
          signature(object = "orderbook"),
          function(object, ...) {
              x <- object@current.ob

              x = x[x[["type"]] == "BID",]

              ## Take care of NAs by using max

              return(max(0, nrow(x)))
          }
          )

## Returns the number of ask orders.

setMethod("ask.orders",
          signature(object = "orderbook"),
          function(object, ...) {
              x <- object@current.ob

              x = x[x[["type"]] == "ASK",]

              return(max(0, nrow(x)))
          }
          )

## Returns the total number of orders.

setMethod("total.orders",
          signature(object = "orderbook"),
          function(object, ...) {

              return(ask.orders(object) + bid.orders(object))

          }
          )

## Returns the midpoint value.

setMethod("mid.point",
          signature(object = "orderbook"),
          function(object, ...) {

              return((best.bid(object)[1] + best.ask(object)[1])/2)

          }
          )

## Returns a data frame with a row for the best ask and a row for the
## best bid.  The columns are price, size, and type.

setMethod("inside.market",
          signature(object = "orderbook"),
          function(object, invis = FALSE, ...) {

              x <- .combine.size(object, .05)

              ## Create ask and bid data frames.

              ask <- x[x[["type"]] == "ASK",]
              bid <- x[x[["type"]] == "BID",]

              ask <- c(price = ask[["price"]][1],
                       size = ask[["size"]][1])

              bid <- c(price = bid[["price"]][nrow(bid)],
                       size = bid[["size"]][nrow(bid)])

              ## If invis is TRUE it won't print but will return the
              ## inside market object

              if(invis == FALSE){
                  cat("Best Bid:          ",
                      .prettify(bid["price"]), "\n")
                  cat("Size:              ",
                      .prettify(bid["size"], "s"), "\n \n")
                  cat("Best Ask:          ",
                      .prettify(ask["price"]), "\n")
                  cat("Size:              ",
                      .prettify(ask["size"], "s"), "\n")
              }

              ## Returns the inside market as an object

              invisible(rbind(ask, bid))
          }
          )

## Returns the spread.

setMethod("spread",
          signature(object = "orderbook"),
          function(object, ...) {

              ask = best.ask(object)
              bid = best.bid(object)

              if(ask[["price"]] == -1 | bid[["price"]] == -1){
                  return(NA)
              } else{
                  return(ask[["price"]] - bid[["price"]])
              }

          }
          )

## Reset to beginning, does not clear trade data or my trades.

setMethod("reset",
          signature(object = "orderbook"),
          function(object){

              ## Clear the hash

              clear(object@ob.data)

              current.ob <- data.frame(numeric(0), numeric(0),
                                       character(0), numeric(0), character(0))

              names(current.ob) <- c("price", "size", "type", "time",
                                     "id")

              object@current.ob <- current.ob
              object@current.time <- 0
              object@file.index <- 1
              object@trade.data <- vector()

              invisible(object)
          }
          )

## Pulls out the orders at the specified price level

setMethod("[",
          signature(x = "orderbook", i = "character"),
          function(x, i){

              current.ob <- x@current.ob
              i = as.numeric(i)

              tmp <- current.ob[current.ob$price == i,]
              rownames(tmp) <- NULL
              return(tmp)
          }
          )

## Creates a copy of the orderbook. This is necessary since hash
## requires copies to be made.

setMethod("copy",
          signature(x = "orderbook"),
          function(x){
              if(length(x@ob.data) > 0)
                  x@ob.data <- copy(x@ob.data)
              else
                  x@ob.data <- hash()

              invisible(x)
          }
          )

## Midpoint Return, automatically finds the midpoint return for the
## selected message row number for a vector of time in seconds,
## e.g. c(5, 10, 60, 120) means find the midpoint return for 5s, 10s,
## 1 min, 2 min after the trade.

setMethod("midpoint.return",
          signature(object = "orderbook"),
          function(object, tradenum, time){

              trade.data <- object@trade.data

              if(isTRUE(object@trader))
                  skip <- 6
              else
                  skip <- 5

              trdprice <- as.numeric(trade.data[4 + (tradenum - 1) * skip])
              trdrow <- as.numeric(trade.data[1 + (tradenum - 1) * skip])

              midpoint.return <- .midpoint.returns(object, trdprice,
                                           trdrow, time)

              midpoint.return <- cbind(midpoint.return)
              rownames(midpoint.return) <- paste(time, "second")

              return(midpoint.return)

          }
          )

## Add an order, you need price, size, type. You should probably
## specify ID but if you don't the orderbook will add it for you
## anyways and automatically make it 1 greater than the current
## greatest ID.

setMethod("add.order",
          signature(object = "orderbook"),
          function(object, price, size, type, time = NULL, id = NULL,
                   status = FALSE){

              ## Some checks to make sure the order is
              ## valid. Specifically, making sure the price and size
              ## are positive and a valid type has been specified.

              stopifnot(price > 0 & size > 0)
              stopifnot(type == "ASK" | type == "BID")

              ## Pull out the order book data frame.

              x <- object@current.ob

              ## If user doesn't specify a time make the new time
              ## 1000ms after the current time

              if(is.null(time)){
                  new.time <- object@current.time + 1000
              } else {
                  new.time <- time
              }

              ## If user doesn't specify an ID make the new id 1
              ## larger than the largest ID.

              if(is.null(id) & nrow(x) != 0){
                  id <- max(as.numeric(x[["id"]])) + 1
              } else if(is.null(id)){
                  id <- 1
              }

              ## Create the new order as a data frame and name it.

              if(isTRUE(object@trader)){
                  new.order <- data.frame(price, size, type, new.time, id, status)
                  names(new.order) <- c("price", "size", "type", "time",
                                        "id", "status")
              } else{
                  new.order <- data.frame(price, size, type, new.time, id)
                  names(new.order) <- c("price", "size", "type", "time",
                                        "id")
              }

              ## Rbind it to current.ob.

              x <- do.call(rbind, list(x, new.order))

              ## Store it into object and return the object.

              object@current.ob <- x
              object@current.time <- new.time

              invisible(object)

          }
          )

## Replace the size of an order.. You need to specify ID and new size.

setMethod("replace.order",
          signature(object = "orderbook"),
          function(object, id, size){

              ## If size is 0 just remove the order.

              if(identical(size, 0)){

                    invisible(remove.order(object, id))

              } else {

                  ## Pull out the order book data frame.

                   x <- object@current.ob

                   ## Make sure the new size isn't greater than the
                   ## current size.

                   tmp.size <- x[["size"]][x[["id"]] == id]
                   if(tmp.size < size){

                       print("Warning size greater than current size")

                   } else {

                       ## Do the replacement.

                       x[["size"]][x[["id"]] == id] <- min(size,
                                    tmp.size)

                       ## Store the new order book data frame then
                       ## return the object.

                       object@current.ob <- x
                       invisible(object)

                   }
               }
          }
          )

## Runs a market order. If there is not enough volume to fill the
## order will be partially filled and cancelled. Might be removed.

setMethod("market.order",
          signature(object = "orderbook"),
          function(object, size, type){

              stopifnot(type == "BUY" | type == "SELL")
              stopifnot(size > 0)

              x <- object@current.ob

              ## Take out ask and bid dataframes.

              ask <- x[x[["type"]] == "ASK",]
              bid <- x[x[["type"]] == "BID",]

              ## If its a buy, then remove orders until you run out
              ## then replace the last one.

              if(type == "BUY" & nrow(ask) > 0){

                  ## Order ask by lowest price, time

                  ask <- ask[order(ask[["price"]], ask[["time"]]),]

                  while(size > 0 & nrow(ask) > 0){

                      ## Decrement size.

                      size = size - ask[["size"]][1]

                      ## Either remove or replace.

                      if(size >= 0){

                          object <- remove.order(object,
                                                 ask[["id"]][1])

                          ask <- ask[-1,]

                      } else if(size < 0){

                          object <- replace.order(object,
                                                  ask[["id"]][1],
                                                  abs(size))
                      }
                  }

                  ## See above

              } else if(type == "SELL" & nrow(bid) > 0){

                  bid <- bid[order(bid[["time"]]),]

                  bid <- bid[order(bid[["price"]], decreasing =
                                   TRUE),]

                  while(size > 0 & nrow(bid) > 0){

                      size <- size - bid[["size"]][1]

                      if(size >= 0){

                          object <- remove.order(object,
                                                 bid[["id"]][1])

                          bid <- bid[-1,]

                      } else if(size < 0){

                          object <- replace.order(object,
                                                  bid[["id"]][1],
                                                  abs(size))
                      }
                  }
              }
              invisible(object)
          }
          )

## Remove an order by ID.

setMethod("remove.order",
          signature(object = "orderbook"),
          function(object, id){

              x <- object@current.ob

              ## Remove from current.ob

              x <- x[x[["id"]] != id,]

              ## Store it back into the object and return the object.

              object@current.ob <- x
              invisible(object)
          }
          )

## View a trade by number

setMethod("view.trade",
          signature(object = "orderbook"),
          function(object, tradenum){

              x <- object@trade.data

              if(isTRUE(object@trader)){
                  skip <- 6
                  names <- c("row", "time", "id", "price", "size", "status")
              }else{
                  skip <- 5
                  names <- c("row", "time", "id", "price", "size")
              }


              start <- (tradenum - 1) * skip + 1
              end <- start + skip - 1

              trade <- x[start:end]
              trade[2] <- .to.time(as.numeric(trade[2]))
              trade <- as.data.frame(trade)

              names(trade) <- paste("trade", tradenum)
              rownames(trade) <- names

              return(trade)
          }
          )
