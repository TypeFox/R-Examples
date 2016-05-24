################################################################################
##
##
## orderbook.plot.R: Helper function for the plot method
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


## Plots the orderbook object at current time. Displays Bids on the left and
## Asks on the right with Price and Size on the Y- and X-axes, respectively.
## Only prices within 10% above and below the midpoint value are shown.

.plot.ob <- function(object, bounds){

    ## Use combine size to find the total size at each price level. This
    ## function returns a data frame. Also get the names for the columns.

    x <- .combine.size(object, bounds)

    ## If there is nothing on the orderbook, stop
    stopifnot(nrow(x)>0)

    ## Maximum size, max/min price and difference between the max
    ## and min price for purposes of drawing the axes.

    max.size <- max(x[["size"]])

    min.price <- min(x[["price"]])
    max.price <- max(x[["price"]])
    midpoint <- mid.point(object)

    bestbid <- best.bid(object)[[1]]
    bestask <- best.ask(object)[[1]]

    ## Creating the x axis values.

    x.at <- pretty(c(0, max.size))
    x.limits <- list(c(x.at[length(x.at)], 0),
                     c(0, x.at[length(x.at)]))

    ## Creating the y axis values and appending the best ask/bestbid to them.

    tmp.at <- c(pretty(c(min.price, max.price), n = 10), bestbid, bestask)
    tmp.at <- sort(tmp.at)

    yask.at <- tmp.at[tmp.at > midpoint]
    ybid.at <- tmp.at[tmp.at < midpoint]

    ## Function for drawing the y-axis.

    new.yscale.components <- function(...) {
        ans <- yscale.components.default(...)
        ans$right <- ans$left

        ans$left$ticks$at <- ybid.at
        ans$left$labels$at <- ybid.at
        ans$left$labels$labels <- formatC(ybid.at, format = "f",
                                          digits = 2)

        ans$right$ticks$at <- yask.at
        ans$right$labels$at <- yask.at
        ans$right$labels$labels <- formatC(yask.at, format = "f",
                                           digits = 2)
        ans
    }

    ## Remove values if its too close to the best bid or best ask.

    space = (max.price - min.price)/20

    if(yask.at[1] + space > yask.at[2])
        yask.at = yask.at[-2]

    if(ybid.at[length(ybid.at)] - space < ybid.at[length(ybid.at) - 1])
        ybid.at = ybid.at[-(length(ybid.at) - 1)]

    ## Ordering the levels so Bid comes before Ask, this allows Bid to be
    ## on the left.

    x[["type"]] <- ordered(x[["type"]], levels = c("BID",
                                                  "ASK"))

    ## Actually plotting it.

    tmp <- xyplot(price~size|type, data = x,

                  ylab = "Price", xlab = "Size (Shares)",

                  main = paste("Order Book",
                  .to.time(object@current.time), sep = " - "),

                  yscale.components = new.yscale.components,

                  scales = list(x = list(relation = "free",
                                limits = x.limits,
                                at = x.at,
                                axs = "i",
                                rot = 45),
                  y = list(alternating = 3)),
                  panel = function(...){
                      panel.xyplot(...)
                      panel.lines(..., type = "H")
                  }
                  )

    ## Return the Trellis object.

    invisible(tmp)

}

## Plot top ask vs top bid, 2nd best ask vs 2nd best bid, etc.

.plot.side.ob <- function(object, n){

    x <- .combine.size(object, 1)
    midpoint = mid.point(object)

    ## Creating the data frame to be plotted.

    ask <- x[x[["type"]] == "ASK",]
    k <- min(n, nrow(ask))

    ask <- ask[1:k,]

    bid <- x[x[["type"]] == "BID",]
    k <- min(n, nrow(bid))

    bid <- bid[(nrow(bid) - k + 1):nrow(bid),]

    x <- rbind(ask, bid)

    y <- data.frame(price = c(seq(ask[1,1], ask[1,1] + (n-1)/100, .01),
                    seq(bid[nrow(bid),1],
                        bid[nrow(bid),1] - (n-1)/100, -.01)),
                    y = c(seq(n, 1, -1)))

    y$price <- round(y$price, 2)

    x <- merge(x, y, all.y = TRUE)

    ## Setting x-axis limits and labels.
    max.size <- ceiling(max(x[["size"]], na.rm = TRUE))

    x.at <- pretty(c(0, max.size))
    x.limits <- c(0, x.at[length(x.at)])

    ## Creating y-axis tick labels.

    yask.at <- rev(x$price[x$price > midpoint])
    ybid.at <- x$price[x$price < midpoint]



   new.yscale.components <- function(...) {
        ans <- yscale.components.default(...)
        ans$right <- ans$left

        ans$left$labels$at <- seq(n)
        ans$left$labels$labels <- formatC(ybid.at, format = "f",
                                          digits = 2)

        ans$right$labels$at <- seq(n)
        ans$right$labels$labels <- formatC(yask.at, format = "f",
                                           digits = 2)

        ans
    }


    new.par.settings <- list(
    layout.widths <- list(left.padding = 2, right.padding = 5))

    tmp <- barchart(y ~ size, data = x, groups = x$type, auto.key = TRUE,
                    ylab = "Bid Price Levels", xlab = "Size (Shares)",
                    main = "Order Book", par.settings = new.par.settings,
                    yscale.components = new.yscale.components,
                    scales = list(x = list(axs = "i",
                                  limits = x.limits,
                                  at = x.at,
                                  tck = c(1, 0)),
                    y = list(alternating = 3))
                    )

    plot(tmp)

    trellis.focus("panel", 1, 1, clip.off = TRUE, highlight = FALSE)
    grid.text("Ask Price Levels", x = 1.14, rot = 90)
    trellis.unfocus()

}

## Same as plot.ob, except # of orders instead of shares at each price
## level.

.plot.orders.ob <-function(object, bounds){

    x <- object@current.ob

    ## We only want data within our bounds

    x <- x[(x[["price"]] < mid.point(object)*(1+bounds)
            & x[["price"]] > mid.point(object)*(1-bounds)),]

    ## Create data.frame with price level, number of orders, and type

    ask <- x[x[["type"]] == "ASK",]
    bid <- x[x[["type"]] == "BID",]

    ask <- data.frame(table(ask[["price"]]))
    bid <- data.frame(table(bid[["price"]]))

    ask <- cbind(ask, rep("ASK", nrow(ask)))
    bid <- cbind(bid, rep("BID", nrow(bid)))

    names(ask) <- c("price", "Orders", "type")
    names(bid) <- names(ask)

    x <- rbind(ask, bid)
    x[["price"]] <- as.numeric(levels(x[["price"]]))

    ## Maximum orders, max/min price. and difference between the max
    ## and min price for purposes of drawing the axes.

    max.orders <- ceiling(max(x[["Orders"]]))

    min.price <- min(x[["price"]])
    max.price <- max(x[["price"]])
    midpoint <- mid.point(object)

    bestbid <- best.bid(object)[[1]]
    bestask <- best.ask(object)[[1]]

    ## Create x axes/limits.

    x.at <- pretty(c(0, max.orders))
    x.limits <- list(c(x.at[length(x.at)], 0),
                     c(0, x.at[length(x.at)]))

    ## Creating the y axis values.

    tmp.at <- c(pretty(c(min.price, max.price), n = 10), bestbid, bestask)
    tmp.at <- sort(tmp.at)

    yask.at <- tmp.at[tmp.at > midpoint]
    ybid.at <- tmp.at[tmp.at < midpoint]

    new.yscale.components <- function(...) {
        ans <- yscale.components.default(...)
        ans$right <- ans$left

        ans$left$ticks$at <- ybid.at
        ans$left$labels$at <- ybid.at
        ans$left$labels$labels <- formatC(ybid.at, format = "f",
                                          digits = 2)

        ans$right$ticks$at <- yask.at
        ans$right$labels$at <- yask.at
        ans$right$labels$labels <- formatC(yask.at, format = "f",
                                           digits = 2)
        ans
    }

    ## Remove values if its too close to the best bid or best ask

    space = (max.price - min.price)/20

    if(yask.at[1] + space > yask.at[2])
        yask.at = yask.at[-2]

    if(ybid.at[length(ybid.at)] - space < ybid.at[length(ybid.at) - 1])
        ybid.at = ybid.at[-(length(ybid.at) - 1)]

    ## Ordering the levels so Bid comes before Ask.

    x[["type"]] <- ordered(x[["type"]],
                                levels = c("BID", "ASK"))

    ## Actually plotting it.

    tmp <- xyplot(x[["price"]]~x[["Orders"]]|x[["type"]], data = x,

                  ylab = "Price", xlab = "Number of Orders",

                  main = paste("Order Book",
                  .to.time(object@current.time), sep = " - "),

                  scales = list(x = list(relation = "free",
                                limits = x.limits,
                                at = x.at,
                                axs = "i"),
                  y = list(alternating = 3)),

                  yscale.components = new.yscale.components,

                  panel = function(...){
                      panel.xyplot(...)
                      panel.lines(..., type = "H")
                  }
                  )

    ## Return the Trellis object

    invisible(tmp)

}



.animate.plot <- function(object, bounds){

    x <- .combine.size(object, 1)
    mid <- mid.point(object)

    ## Find the min ask and max bid price for this current.ob

    ask <- x[x[["type"]] == "ASK",]
    ask <- ask[ask$price < min(ask$price) + bounds - .001,]

    bid <- x[x[["type"]] == "BID",]
    bid <- bid[bid$price > max(bid$price) - bounds + .001,]

    ## Set y limits

    y.limits <- c(min(bid$price), max(ask$price))

    ## Find the max size for this current.ob

    max.size <- max(x$size[x$price <= y.limits[2] & x$price >=
                               y.limits[1]])

    x.at <- pretty(c(0, max.size))
    x.limits <- list(c(x.at[length(x.at)], 0),
                     c(0, x.at[length(x.at)]))

    x <- object@current.ob

    ask <- x[x[["type"]] == "ASK",]
    ask <- ask[ask$price < y.limits[2] + .001,]

    bid <- x[x[["type"]] == "BID",]
    bid <- bid[bid$price > y.limits[1] - .001,]

    x <- do.call(rbind, list(ask, bid))
    x <- x[order(x$price),]

    price <- round(seq(y.limits[1], y.limits[2], .01), 2)

    price <- cbind(price)

    x <- merge(x, price, all.y = TRUE)

    x$price <- as.ordered(x$price)
    x$time <- as.ordered(x$time)

    ## Ordering the levels so Bid comes before Ask.

    x[["type"]] <- ordered(x[["type"]],
                           levels = c("BID", "ASK"))

    ## Creating the y-axis

    ymin = y.limits[1] - .01
    ybid.at <- (100 * (min(bid$price) - ymin)):(100 * (max(bid$price) - ymin))
    ybid.at <- round(ybid.at)
    yask.at <- (100 * (min(ask$price) - ymin)):(100 * (max(ask$price) - ymin))
    yask.at <- round(yask.at)

    new.yscale.components <- function(...) {
        ans <- yscale.components.default(...)
        ans$right <- ans$left

        ans$left$ticks$at <- ybid.at
        ans$left$labels$at <- ybid.at
        ans$left$labels$labels <- formatC(price[ybid.at], format =
                                          "f", digits = 2)

        ans$right$ticks$at <- yask.at
        ans$right$labels$at <- yask.at
        ans$right$labels$labels <- formatC(price[yask.at], format =
                                           "f", digits = 2)

        ans
    }

    if(isTRUE(object@trader)){
        groups <- interaction(x$status, x$time)
        colors <- c("black", "gray")
    }else{
        groups <- x$time
        colors <- "gray"
    }

    ## Actually plotting it

    tmp <- barchart(price ~ size | type, data = x,

                    ylab = "Price", xlab = "Size (Shares)",
                    groups = groups,

                    main = paste("Order Book",
                    .to.time(object@current.time), sep = " -- "),

                    stack = TRUE,

                    col = colors,


                    scale = list(x = list(relation = "free", at = x.at,
                                 limits = x.limits, axs = "i", rot = 45),
                    y = list(alternating = 3)),
                    yscale.components = new.yscale.components
                    )

    invisible(tmp)
}


.supply.demand.plot <- function(object, bounds){

    x <- .combine.size(object, bounds)
    ask <- x[x$type == "ASK",]
    bid <- x[x$type == "BID",]
    bid <- bid[order(bid$price, decreasing = TRUE),]

    ask$price <- ask$price - mid.point(object)
    bid$price <- bid$price - mid.point(object)
    ask$price <- ask$price/max(ask$price)
    bid$price <- bid$price/max(abs(bid$price))

    ask$size <- cumsum(ask$size)/sum(ask$size)
    bid$size <- cumsum(bid$size)/sum(bid$size)
    x <- rbind(ask, bid)
    x <- rbind(x, c(ask$price[1], 0, "ASK"), c(bid$price[1], 0, "BID"))
    x$size <- as.numeric(x$size)
    x$price <- as.numeric(x$price)

    x.limits = c(0, 1.2)
    x.at = c(0, .2, .4, .6, .8, 1, 1.2)
    x.labels = c(0, 20, 40, 60, 80, 100, 120)

    y.limits <- c(-1.5, 1.5)
    y.at = c(-1.5, -1, -.5, 0, .5, 1, 1.5)
    y.labels = c(-150, -100, -50, 0, 50, 100, 150)

    tmp <- xyplot(x$price ~ x$size, data = x, groups = x$type, type =
                  "S", ylab = "Price (%)",
                  xlab = "Size (%)", main = "Supply and Demand", sub =
                  .to.time(object@current.time),
                  scales = list(

                  x = list(limits = x.limits, axs = "i", at = x.at, labels =
                  x.labels, tck = c(1, 0)),

                  y = list(limits = y.limits, at = y.at, labels = y.labels, tck =
                  c(1, 0))),

                  panel = function(...){
                      panel.xyplot(...)
                      panel.abline(h = 0)
                  }
                  )

    invisible(tmp)

}

