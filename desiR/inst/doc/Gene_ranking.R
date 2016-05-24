## ------------------------------------------------------------------------
library(desiR)

## ----fig.show='hold',  fig.height=4, fig.width=6-------------------------
# generate values to convert to desirabilities
x <- seq(0,1, length.out=500)

# view mappings of x to different desirabilities
par(mfrow=c(2,3),
    mar=c(3,4.5,3,1),
    las=1)

# "high is good"
d1 <- d.high(x, cut1 = 0.2, cut2 = 0.8)
plot(d1 ~ x, ylab="Desirability", ylim=c(0,1), type="l",
     main="High")

# "low is good", with non-default scale and min.des arguments
d2 <- d.low(x, cut1 = 0.2, cut2 = 0.8, scale = 0.5, des.min = 0.1)
plot(d2 ~ x, ylab="Desirability", ylim=c(0,1), type="l",
     main="Low")

# "high is good" with a four parameter logistic function
d3 <- d.4pl(x, hill = 5, inflec = 0.5)
plot(d3 ~ x, ylab="Desirability", ylim=c(0,1), type="l",
     main="Logistic")

# "central or middle values are good"
d4 <- d.central(x, cut1 = 0.1, cut2 = 0.2, cut3 = 0.6, cut4 = 0.8)
plot(d4 ~ x, ylab="Desirability", ylim=c(0,1), type="l",
     main="Central")

# "extreme values are good"
d5 <- d.ends(x, cut1 = 0.2, cut2 = 0.3, cut3 = 0.7, cut4 = 0.8)
plot(d5 ~ x, ylab="Desirability", ylim=c(0,1), type="l",
     main="Ends")

# categorical: C is most desirable, followed by B
barplot(c(0.1, 0.5, 1), ylim=c(0,1.1),ylab="Desirability",
        main="Categorical", names.arg=c("A","B","C")); box()

## ------------------------------------------------------------------------
# load the data (provided with the package)
data(farmer2005)

# look at first few rows
head(farmer2005)

## ---- fig.show='hold', fig.height=8, fig.width=5-------------------------
par(mfrow=c(3,2),
    las=1)
hist(farmer2005$AveExpr, breaks=30, col="grey", border="white", main="",
     xlab="Mean expression")
des.line(farmer2005$AveExpr, "d.high", des.args=c(cut1=5.75, cut2=7, scale=0.5))

hist(farmer2005$SD, breaks=30, col="grey", border="white", main="",
     xlab="Within group SD")
des.line(farmer2005$SD, "d.high", des.args=c(cut1=0.1, cut2=0.3, scale=0.5))

hist(farmer2005$P.Value, breaks=50, col="grey", border="white", main="",
     xlab="p-value")
des.line(farmer2005$P.Value, "d.low", des.args=c(cut1=0.0001, cut2=0.1, scale=2))

hist(farmer2005$logFC, breaks=50, col="grey", border="white", main="",
     xlab="Log FC")
des.line(farmer2005$logFC, "d.ends", des.args=c(cut1=log2(1/3), cut2=log2(1/1.5),
         cut3=log2(1.5), cut4=log2(3), scale=0.5))

hist(abs(farmer2005$PCNA.cor), breaks=50, col="grey", border="white", main="",
     xlab="Absolute correlation with PCNA")
des.line(farmer2005$P.Value, "d.low",  des.args=c(cut1=0.15, cut2=0.5,
         scale=1.5, des.min = 0.05))

## ------------------------------------------------------------------------
farmer2005$d1 <- d.high(farmer2005$AveExpr, cut1=5.75, cut2=7, scale=0.5)
farmer2005$d2 <- d.high(farmer2005$SD, cut1=0.1, cut2=0.3, scale=0.5)
farmer2005$d3 <- d.low(farmer2005$P.Value, cut1=0.0001, cut2=0.1, scale=2)
farmer2005$d4 <- d.ends(farmer2005$logFC, cut1=log2(1/3), cut2=log2(1/1.5),
         cut3=log2(1.5), cut4=log2(3), scale=0.5)
farmer2005$d5 <- d.low(abs(farmer2005$PCNA.cor), cut1=0.15, cut2=0.5,
            scale=1.5, des.min = 0.05)

## ------------------------------------------------------------------------
# order of weights needs to match the order of desirabilities
farmer2005$D <- d.overall(farmer2005$d1, farmer2005$d2, farmer2005$d3,
                          farmer2005$d4, farmer2005$d5,
                weights=c(0.1, 0.1, 1, 1, 0.5))

## ----fig.height=4, fig.width=5-------------------------------------------
par(las=1)
plot(rev(sort(farmer2005$D)), type="l", xlab="Rank", ylab="Overall Desirability")

## ------------------------------------------------------------------------
f2 <- farmer2005[order(farmer2005$D, decreasing=TRUE), ]
head(f2, 10)

## ---- fig.height=4, fig.width=4------------------------------------------
# volcano plot
par(las=1)
plot(-log10(P.Value) ~ logFC, data=farmer2005, col="darkgrey", cex=0.75)

# overplot the top 8 genes
points(-log10(P.Value) ~ logFC, pch=16, 
       data=farmer2005[farmer2005$D > 0.8, ])

# horizontal reference line at p = 0.05.
abline(h=-log10(0.05), lty=2)

