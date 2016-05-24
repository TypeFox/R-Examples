# code for Figure 7.5


library(MASS)

# make dataset

dfr1 = data.frame(mvrnorm(n=100, rep(6, 2), matrix(c(10,3,3,2),2,2)))
colnames(dfr1) = c("x", "y")
dfr2 = dfr1
dfr2 = dfr2[order(dfr2$x),]
delta = 2
dfr2[1,2] = dfr2[1,2]+delta
dfr2[nrow(dfr2),2] = dfr2[nrow(dfr2),2]-delta
dfr2 = dfr2[order(as.numeric(rownames(dfr2))),]
dfr1$x1 = dfr2$x
dfr1$y1 = dfr2$y
dfr2 = dfr1
dfr2$x = dfr2$x - mean(dfr2$x)
dfr2$x1 = dfr2$x1 - mean(dfr2$x1)
dfr2$which = "centered"
dfr1$which = "uncentered"
dfr = rbind(dfr1, dfr2)
dfr$which = relevel(factor(dfr$which), "uncentered")

# plot with lattice

xyplot(y~x|which, data=dfr, x1=dfr$x1, y1=dfr$y1,
  panel = function(x, y, x1, y1, subscripts, ...) {
    x1 <- x1[subscripts]
    y1 <- y1[subscripts]
    panel.grid(h = -1, v = -1)
    panel.points(x, y)
    panel.abline(lm(y~x))
    panel.abline(lm(y1~x1), lty=2, col.line="black")
    panel.abline(v=0, col.line="darkgrey")
  }
)

# or simply with plot():

par(mfrow=c(1,2))
plot(dfr[dfr$which=="uncentered",]$x, dfr[dfr$which=="uncentered",]$y)
abline(lm(y~x, data= dfr, subset=which=="uncentered"))
abline(lm(y1~x1, data= dfr, subset=which=="uncentered"))
abline(v=0)

plot(dfr[dfr$which=="centered",]$x, dfr[dfr$which=="centered",]$y)
abline(lm(y~x, data= dfr, subset=which=="centered"))
abline(lm(y1~x1, data= dfr, subset=which=="centered"))
abline(v=0)
par(mfrow=c(1,1))

