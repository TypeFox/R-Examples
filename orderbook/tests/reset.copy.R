library(orderbook)
load("read.orders.time.test.RData")

filename <- system.file("extdata", "sample.txt",
                        package = "orderbook")

ob <- orderbook(file = filename)
ob <- read.orders(ob, 5000)

## Test copy

ob.test <- copy(ob)

stopifnot(isTRUE(identical(ob@current.ob, ob.test@current.ob)))

## Test reset

ob <- reset(ob)

stopifnot(isTRUE(identical(0, ob@current.time)))
stopifnot(isTRUE(identical(1, ob@file.index)))
stopifnot(isTRUE(identical(as.integer(0), nrow(ob@current.ob))))
stopifnot(isTRUE(identical(as.integer(0), length(ob@ob.data))))
