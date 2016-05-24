library(orderbook)

filename <- system.file("extdata", "sample.txt",
                        package = "orderbook")

ob <- orderbook(file = filename)
ob <- read.orders(ob, 5000)
ob <- next.trade(ob)

## Next trade

stopifnot(isTRUE(identical(5045, ob@file.index)))

## Trade after

ob <- next.trade(ob)

stopifnot(isTRUE(identical(5047, ob@file.index)))

## Trade before

ob <- previous.trade(ob)

stopifnot(isTRUE(identical(5045, ob@file.index)))

## Trade before

ob <- previous.trade(ob)

stopifnot(isTRUE(identical(4994, ob@file.index)))
