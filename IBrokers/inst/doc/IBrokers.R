### R code from vignette source 'IBrokers.Rnw'

###################################################
### code chunk number 1: twsconnect (eval = FALSE)
###################################################
## library(IBrokers)
## tws <- twsConnect()
## tws
## reqCurrentTime(tws)
## serverVersion(tws)
## twsDisconnect(tws)


###################################################
### code chunk number 2: twsconnect
###################################################
if(Sys.getenv("IBCON")=="1") {
library(IBrokers)
tws <- twsConnect(999, port=Sys.getenv("IBPORT"))
tws
reqCurrentTime(tws)
serverVersion(tws)
twsDisconnect(tws)
}


###################################################
### code chunk number 3: reqmktdata (eval = FALSE)
###################################################
## tws <- twsConnect()
## twsFuture("YM","ECBOT","200809")
## reqContractDetails(tws, twsEquity("QQQQ"))


###################################################
### code chunk number 4: reqmktdata
###################################################
if(Sys.getenv("IBCON")=="1") {
tws <- twsConnect(999, port=Sys.getenv("IBPORT"))
twsFuture("YM","ECBOT","200809")
reqContractDetails(tws, twsEquity("QQQQ"))
}


###################################################
### code chunk number 5: reqmktdata2 (eval = FALSE)
###################################################
## reqMktData(tws, twsEquity("QQQQ"))


###################################################
### code chunk number 6: reqmktdata3
###################################################
if(Sys.getenv("IBCON")=="1") {
IBrokers:::.reqMktData.vignette(tws, twsEquity("QQQQ"))
twsDisconnect(tws)
}


###################################################
### code chunk number 7: playbackmktdata (eval = FALSE)
###################################################
## reqMktData(tws, twsEquity("SBUX"), CALLBACK=NULL, file="SBUX.dat")
## twsp <- twsConnect(filename="SBUX.dat")


###################################################
### code chunk number 8: playbackmktdata
###################################################
if(Sys.getenv("IBCON")=="1") {
tws <- twsConnect(999, port=Sys.getenv("IBPORT"))
IBrokers:::.reqMktData.vignette(tws, twsEquity("SBUX"), CALLBACK=NULL, file="SBUX.dat")
twsDisconnect(tws)
twsp <- twsConnect(filename="SBUX.dat")
}


###################################################
### code chunk number 9: playback (eval = FALSE)
###################################################
## reqMktData(twsp)
## reqMktData(twsp, playback=0)


###################################################
### code chunk number 10: playback
###################################################
if(Sys.getenv("IBCON")=="1") {
reqMktData(twsp)
reqMktData(twsp, playback=0)
}


###################################################
### code chunk number 11: closeplayback
###################################################
if(Sys.getenv("IBCON")=="1") {
twsDisconnect(twsp)
}


