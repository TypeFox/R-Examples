library(orderbook)

filename <- system.file("extdata", "sample.txt",
                        package = "orderbook")

ob <- orderbook(file = filename)

## Add order Bid

ob <- add.order(ob, price = 123, size = 123, type = "BID", time = 1, id = 1)

## Add order Ask

ob <- add.order(ob, price = 124, size = 124, type = "ASK")


stopifnot(isTRUE(1001 %in% ob@current.ob$time))
stopifnot(isTRUE(identical(123.5, mid.point(ob)[[1]])))

## Add more orders

ob <- add.order(ob, price = 125, size = 124, type = "ASK")
ob <- add.order(ob, price = 126, size = 124, type = "ASK")

## Remove an order

ob <- remove.order(ob, 3)

stopifnot(!isTRUE(3 %in% ob@current.ob$id))

## Replace the size of the bid

ob <- replace.order(ob, 1, 100)

stopifnot(isTRUE(100 %in% ob@current.ob$size))

## Do a market order

ob <- market.order(ob, 125, "BUY")

stopifnot(!isTRUE(2 %in% ob@current.ob$id))
stopifnot(isTRUE(123 %in% ob@current.ob$size))
