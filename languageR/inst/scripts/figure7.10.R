# figure 7.10

data(beginningReaders)
beginningReaders$OrthLength = scale(beginningReaders$OrthLength, scale=FALSE)
beginningReaders$LogFrequency = scale(beginningReaders$LogFrequency,scale=FALSE)


sm1 = lmer(LogRT ~ PC1+PC2+PC3 + ReadingScore +
OrthLength + I(OrthLength^2) + LogFrequency + LogFamilySize +
(1|Word) + (1|Subject) + (0+LogFrequency|Subject) + (0+OrthLength|Subject), 
data=beginningReaders)

df <- coef(lmList(LogRT ~ PC1+PC2+PC3 + ReadingScore +
OrthLength + I(OrthLength^2) + LogFrequency + LogFamilySize | Subject, 
data=beginningReaders))[,c(1,8)]

cc1 <- data.frame(A = coef(sm1)[[2]][,1], 
                  B = coef(sm1)[[3]][,8])
names(cc1) <- c("A", "B")
df <- cbind(df, cc1)
with(df,
  print(xyplot(`(Intercept)` ~ LogFrequency, aspect = 1,
    x1 = B, y1 = A, 
    panel = function(x, y, x1, y1, subscripts, ...) {
    panel.grid(h = -1, v = -1)
      x1 <- x1[subscripts]
      y1 <- y1[subscripts]
      panel.abline(v=summary(sm1)@coefs["LogFrequency",1],col="darkgrey",lwd=2)
      panel.abline(h=summary(sm1)@coefs["(Intercept)",1],col="darkgrey",lwd=2)
      panel.arrows(x, y, x1, y1, type = "closed", length = 0.1,
        angle = 15, ...)
      panel.points(x, y,
        col = trellis.par.get("superpose.symbol")$col[2])
      panel.points(x1, y1, pch = 3)
    },
    key = list(space = "top", columns = 2,
    text = list(c("Mixed model", "Within-group")),
    points = list(col = trellis.par.get("superpose.symbol")$col[1:2],
                  pch = c(3,1)
    ))
)))

