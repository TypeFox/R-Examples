library(orderbook)
load("read.orders.time.test.RData")

filename <- system.file("extdata", "sample.txt",
                        package = "orderbook")

ob <- orderbook(file = filename)
ob <- read.orders(ob, 5000)

## Best bid

stopifnot(isTRUE(identical(11.4, best.bid(ob)[[1]])))

## Best ask

stopifnot(isTRUE(identical(11.41, best.ask(ob)[[1]])))

## Bid Price Levels

stopifnot(isTRUE(identical(as.integer(181), bid.price.levels(ob))))

## Ask Price Levels

stopifnot(isTRUE(identical(as.integer(534), ask.price.levels(ob))))

## Total Price Levels

stopifnot(isTRUE(identical(as.integer(715), total.price.levels(ob))))

## Bid orders

stopifnot(isTRUE(identical(615, bid.orders(ob))))

## Ask orders

stopifnot(isTRUE(identical(1834, ask.orders(ob))))

## Total orders

stopifnot(isTRUE(identical(2449, total.orders(ob))))

## Midpoint

stopifnot(isTRUE(identical(11.405, round(mid.point(ob)[[1]], 3))))

## Inside.market

test <- inside.market(ob, invis = TRUE)

stopifnot(isTRUE(identical(11.41, test[[1]])))
stopifnot(isTRUE(identical(11.4, test[[2]])))
stopifnot(isTRUE(identical(300, test[[3]])))
stopifnot(isTRUE(identical(7400, test[[4]])))

## Spread

stopifnot(isTRUE(identical(0.01, round(spread(ob), 2))))

## "["

ids = c("6487521", "6552654", "6598723")

stopifnot(isTRUE(all(ids %in% ob["11.41"]$id)))

## Midpoint Returns test

midpoints <- midpoint.return(ob, 5, c(5, 10))

stopifnot(isTRUE(identical(round(midpoints[[1]], 2), 0.02)))
stopifnot(isTRUE(identical(round(midpoints[[2]], 2), 0.02)))
