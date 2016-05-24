setOldClass("mts")
setOldClass("ts")
setOldClass("xts")
setOldClass("zoo")
##
## Class definition for portfolio solution
##
setClass("PortSol", representation(weights = "numeric", opt = "list", type = "character", call = "call"))
##
## Class definition for portfolio CDaR
##
setClass("PortCdd", representation(CDaR = "numeric", thresh = "numeric", DrawDown = "timeSeries"), contains = "PortSol")
##
## Class definition for portfolio aveDD
##
setClass("PortAdd", representation(AveDD = "numeric", DrawDown = "timeSeries"), contains = "PortSol")
##
## Class definition for portfolio maxDD
##
setClass("PortMdd", representation(MaxDD = "numeric", DrawDown = "timeSeries"), contains = "PortSol")
##
## Defining Class Union for draw down portfolios
## 
setClassUnion("PortDD", c("PortAdd", "PortCdd", "PortMdd"))

