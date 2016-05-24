library(QCA3)
data(CarenPanofsky)
tqca.tt <- cs_truthTable(CarenPanofsky,'recognition',names(CarenPanofsky)[1:5])
tqca.ans <- reduce(tqca.tt) 
QCA3:::prettyPI(tqca.ans)

data(McCammonVanDyke)
workdat <- McCammonVanDyke
workdat[workdat==-9] <- 0
fig13.2 <- reduce(workdat,"coalition",c("ideology","threats","opportunity","ties","resources"))
QCA3:::prettyPI(fig13.2)
## result in figure 13.2
     
workdat <- McCammonVanDyke
idx <- apply(workdat, 1, function(x) any(x==-9))
ans <- reduce(workdat[!idx,],"coalition",c("ideology","threats","opportunity","ties","resources"))
fig13.3 <- constrReduce(ans,include=workdat[idx,1:5])
QCA3:::prettyPI(fig13.3)
