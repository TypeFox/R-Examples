### R code from vignette source 'portfolioSim.Rnw'

###################################################
### code chunk number 1: portfolioSim.Rnw:24-28
###################################################
## Sets display options
options(width = 75, digits = 2, scipen = 5)

library(portfolioSim)


###################################################
### code chunk number 2: portfolioSim.Rnw:124-125
###################################################
data.frame(id = c("A", "B", "C"), shares = c(10, 20, -10))


###################################################
### code chunk number 3: portfolioSim.Rnw:140-141
###################################################
data.frame(id = c("B","C","D","E"), side=c("S","C","B","X"), shares=c(20,5,10,20))


###################################################
### code chunk number 4: portfolioSim.Rnw:147-148
###################################################
data.frame(id=c("A","C","D","E"), shares=c(10,-5,10,-20))


###################################################
### code chunk number 5: portfolioSim.Rnw:257-298
###################################################

library(portfolioSim)

sdi <- new("sdiDf", data = data.frame(period = c(rep(1, times=4),
                                                 rep(2, times=4),
                                                 rep(3, times=4),
                                                 rep(4, times=4)),
                                      id = rep(c("A","B","C","D"), times=4)))
sdi@data$id <- as.character(sdi@data$id)
sdi@data$ret <- 0
sdi@data$start.price <- 10
sdi@data$end.price <- 10
sdi@data$volume <- 100
sdi@data$universe <- TRUE
sdi@data$universe[15] <- FALSE

trades.1 <- new("trades", trades = data.frame(id = c("A","B"),
                                              side = c("B","X"),
                                              shares = c(10,10)))
trades.2 <- new("trades", trades = data.frame(id = c("B","C","A","A"),
                                              side = c("C","B","X","S"),
                                              shares = c(5,5,10,10)))
trades.3 <- new("trades", trades = data.frame(id = c("B","D"),
                                              side = c("C","B"),
                                              shares = c(5,20)))
trades.4 <- new("trades", trades = data.frame(id = c("D"),
                                              side = c("B"),
                                              shares = c(5)))

sti <- new("stiPresetTrades", periods = c(1,2,3,4),
                              sim.trades = list(trades.1, trades.2,
                                                trades.3, trades.4))

periods <- data.frame(period = c(1,2,3,4), start = c(1,2,3,4), end = c(1,2,3,4))

ps <- new("portfolioSim", periods = periods, data.interface = sdi,
          trades.interface = sti, out.loc = "out_dir_1",
          out.type = c("portfolio","trades","basic"))

sim.portfolios <- runSim(ps, verbose = FALSE)



###################################################
### code chunk number 6: portfolioSim.Rnw:312-313
###################################################
getSimData(ps@data.interface, period = 1, verbose = FALSE)@data


###################################################
### code chunk number 7: portfolioSim.Rnw:320-322
###################################################
getSimTrades(ps@trades.interface, period = 1, holdings = new("portfolio"),
             sim.data = new("simData"), verbose = FALSE)@trades@trades


###################################################
### code chunk number 8: portfolioSim.Rnw:331-332
###################################################
sim.portfolios@data[[1]]@end.data@holdings@shares


###################################################
### code chunk number 9: portfolioSim.Rnw:340-342
###################################################
getSimTrades(ps@trades.interface, period = 2, holdings = new("portfolio"),
             sim.data = new("simData"), verbose = FALSE)@trades@trades


###################################################
### code chunk number 10: portfolioSim.Rnw:352-353
###################################################
sim.portfolios@data[[2]]@end.data@holdings@shares


###################################################
### code chunk number 11: portfolioSim.Rnw:359-361
###################################################
getSimTrades(ps@trades.interface, period = 3, holdings = new("portfolio"),
             sim.data = new("simData"), verbose = FALSE)@trades@trades


###################################################
### code chunk number 12: portfolioSim.Rnw:370-371
###################################################
sim.portfolios@data[[3]]@end.data@holdings@shares


###################################################
### code chunk number 13: portfolioSim.Rnw:383-385
###################################################
getSimTrades(ps@trades.interface, period = 4, holdings = new("portfolio"),
             sim.data = new("simData"), verbose = FALSE)@trades@trades


###################################################
### code chunk number 14: portfolioSim.Rnw:391-392
###################################################
getSimData(ps@data.interface, period = 4, verbose = FALSE)@data


###################################################
### code chunk number 15: portfolioSim.Rnw:401-402
###################################################
sim.portfolios@data[[4]]@start.data@holdings@shares


###################################################
### code chunk number 16: portfolioSim.Rnw:408-409
###################################################
sim.portfolios@data[[4]]@end.data@holdings@shares


###################################################
### code chunk number 17: portfolioSim.Rnw:455-458
###################################################
data.frame(period = 1:4,
start = as.POSIXct(c("2006-01-01 09:30", "2006-04-01 09:30", "2006-07-01 09:30", "2006-10-01 09:30")),
end = as.POSIXct(c("2006-03-31 16:00", "2006-06-30 16:00", "2006-09-30 16:00", "2006-12-31 16:00")))


###################################################
### code chunk number 18: portfolioSim.Rnw:860-862
###################################################
data(starmine.sim)
names(starmine.sim)


###################################################
### code chunk number 19: portfolioSim.Rnw:869-875
###################################################
periods <- data.frame(period = sort(starmine.sim$date[!duplicated(starmine.sim$date)]))
## periods <- data.frame(period = unique(starmine.sim$date))
periods$start <- periods$period
periods$end <- c(periods$start[-1], as.Date("1995-12-31"))

periods


###################################################
### code chunk number 20: portfolioSim.Rnw:904-905
###################################################
starmine.sim$period <- starmine.sim$date


###################################################
### code chunk number 21: portfolioSim.Rnw:918-920
###################################################
starmine.sim$start.price <- starmine.sim$prior.close.usd
starmine.sim$end.price <- starmine.sim$price.usd


###################################################
### code chunk number 22: portfolioSim.Rnw:940-941
###################################################
starmine.sim$ret <- starmine.sim$ret.1m


###################################################
### code chunk number 23: portfolioSim.Rnw:952-953
###################################################
starmine.sim$universe <- TRUE


###################################################
### code chunk number 24: portfolioSim.Rnw:960-961
###################################################
data.interface <- new("sdiDf", data = starmine.sim)


###################################################
### code chunk number 25: portfolioSim.Rnw:986-988
###################################################
trades.interface <- new("stiFromSignal", in.var = "smi",
size = "decile", equity = 1000000, rebal.on = periods$period)


###################################################
### code chunk number 26: portfolioSim.Rnw:1005-1008
###################################################
ps <- new("portfolioSim", periods = periods, freq = 12,
data.interface = data.interface, trades.interface = trades.interface,
fill.volume.pct = Inf, out.loc = "out_dir_2", out.type = "lean")


###################################################
### code chunk number 27: portfolioSim.Rnw:1026-1027
###################################################
result <- runSim(ps, verbose = FALSE)


###################################################
### code chunk number 28: portfolioSim.Rnw:1047-1048
###################################################
summary(result)


###################################################
### code chunk number 29: portfolioSim.Rnw:1068-1070
###################################################
ps@out.type <- "detail"
result <- runSim(ps, verbose = FALSE)


###################################################
### code chunk number 30: portfolioSim.Rnw:1073-1074
###################################################
summary(result)


###################################################
### code chunk number 31: portfolioSim.Rnw:1081-1082
###################################################
result <- loadIn(new("simResult"), in.loc = "out_dir_2")


###################################################
### code chunk number 32: portfolioSim.Rnw:1156-1158
###################################################
unlink("out_dir_1", recursive = TRUE)
unlink("out_dir_2", recursive = TRUE)


