### R code from vignette source 'orderbook.Rnw'

###################################################
### code chunk number 1: orderbook.Rnw:128-134
###################################################
library(orderbook)
filename <- system.file("extdata",
                        "sample.txt",
                        package = "orderbook")
ob <- orderbook(file = filename)
ob <- read.orders(ob, 10000)


###################################################
### code chunk number 2: orderbook.Rnw:136-137
###################################################
ob


###################################################
### code chunk number 3: orderbook.Rnw:147-148
###################################################
summary(ob)


###################################################
### code chunk number 4: orderbook.Rnw:159-160
###################################################
display(ob)


###################################################
### code chunk number 5: orderbook.Rnw:169-170
###################################################
plot(ob)


###################################################
### code chunk number 6: orderbook.Rnw:182-183
###################################################
ob["11.01"]


###################################################
### code chunk number 7: orderbook.Rnw:197-198
###################################################
nrow(ob["11.00"])


###################################################
### code chunk number 8: orderbook.Rnw:207-208
###################################################
plot(ob, bounds = 0.033, type = 'o')


###################################################
### code chunk number 9: orderbook.Rnw:225-226
###################################################
plot(ob, bounds = 0.01, type = "sd")


###################################################
### code chunk number 10: orderbook.Rnw:237-238
###################################################
ob <- read.time(ob, "9:30:00")


###################################################
### code chunk number 11: orderbook.Rnw:246-248
###################################################
ob <- read.orders(ob, n = -50)
ob


###################################################
### code chunk number 12: orderbook.Rnw:290-297
###################################################
filename <- system.file("extdata",
                        "tradersample.txt",
                        package = "orderbook")
ob <- orderbook(file = filename)
ob <- read.time(ob, "9:30:05")
ob <- next.trade(ob)
ob


###################################################
### code chunk number 13: orderbook.Rnw:306-308
###################################################
view.trade(ob, tradenum = 584)
mid.point(ob)


###################################################
### code chunk number 14: orderbook.Rnw:316-317
###################################################
midpoint.return(ob, tradenum = 584, time = 10)


###################################################
### code chunk number 15: orderbook.Rnw:330-332
###################################################
ob <- read.time(ob, "9:30:15")
plot(ob, type = "t", bounds = 0.02)


###################################################
### code chunk number 16: orderbook.Rnw:354-358 (eval = FALSE)
###################################################
## ob <- add.order(ob, 11.20, 300, "ASK")
## ob <- remove.order(ob, 1231883)
## ob <- replace.order(ob, 1231883, 150)
## ob <- market.order(ob, 200, "BUY")


###################################################
### code chunk number 17: orderbook.Rnw:376-469
###################################################
simulate <- function(ob, n=1000,
                     action.prob = c(cancel=0.5, market=0.2, limit=0.30),
                     order.type  = 0.5,
                     alpha = 0.3,
                     in.spread.prob = 0.35,
                     ...){

    x = ob@current.ob

    tmp.midpoint = mid.point(ob)
    tmp.bestask = best.ask(ob)
    tmp.bestbid = best.bid(ob)

    for(i in 1:n) {

        x = ob@current.ob
        if(mid.point(ob) == 0){
            current.price = tmp.midpoint
        } else {
            current.price <- mid.point(ob)
        }

        isbuy = runif(1) < order.type

        if(total.orders(ob) < 250){
            action.prob[1] <- 0
            action.prob[4] = 1 - sum(action.prob[1:3])
        } else {
            action.prob[1] = 0.5
            action.prob[4] = 1 - sum(action.prob[1:3])
        }


        action <- sample(c("Cancel", "Market", "Limit", ""),
                         size=1, prob=action.prob)

        if (action == "Cancel") {
              ## pick an existing ID and cancel the order

            ob <- remove.order(ob, sample(x[["id"]], size = 1))

        }
        else if (action == "Market") {

            ## set a new price/ or tick
            if(isbuy) {
                ob <- market.order(ob, type="BUY",
                                   size = best.ask(ob)[2] )

            } else {
                  ob <- market.order(ob, type="SELL",
                                     size = best.bid(ob)[2] )
              }

        }
        if (action == "Limit") {

            if(spread(ob) <= 0.01){
                spread.diff = 0
            } else {
                spread.diff = round(runif(1, 0, spread(ob)), 2)
            }

            out.diff = round((mid.point(ob)*.1)*runif(1)^1/(1 + alpha), 2)

            in.spread = runif(1) < in.spread.prob
            size = round(exp(rnorm(1, mean = 4.5, sd = .8)))

              if(isbuy & in.spread){

                  ob <- add.order(ob, price= max(0,
                                      best.bid(ob)[1] + spread.diff),
                                  size, type="BID")
              } else if(isbuy & !in.spread){
                  ob <- add.order(ob, price = max(0,
                                      best.bid(ob)[1] - out.diff),
                                  size, type = "BID")
              } else if(!isbuy & in.spread){

                  ob <- add.order(ob, price= max(0,
                                      best.ask(ob)[1] - spread.diff),
                                  size, type="ASK")
              } else if (!isbuy & !in.spread){
                  ob <- add.order(ob, price = max(0,
                                      best.ask(ob)[1] + out.diff),
                                  size, type = "ASK")
              }
        }
    }


    invisible(ob)
}


###################################################
### code chunk number 18: orderbook.Rnw:472-473
###################################################
ob <- simulate(ob)


###################################################
### code chunk number 19: orderbook.Rnw:478-479
###################################################
plot(ob)


