library("R2BayesX")

## load zambia data and map
data("ZambiaNutrition")
data("ZambiaBnd")

## estimate model
zm1 <- bayesx(stunting ~ memployment + meducation + urban + gender + 
  sx(mbmi) + sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) +
  sx(district, bs = "re"), iter = 12000, burnin = 2000, step = 10,
  data = ZambiaNutrition)

## summary statistics
summary(zm1)

## visualising estimation results
plot(zm1, term = "sx(mbmi)")
plot(zm1, term = "sx(agechild)")

plot(zm1, term = "sx(district)")
plot(zm1, term = c("sx(district)"), 
  map = ZambiaBnd, pos = "topleft")

## customizing graphics
plot(zm1, term = "sx(mbmi)", resid = TRUE, cex = 0.1)
plot(zm1, term = "sx(agechild)", resid = TRUE, cex = 0.1)

plot(zm1, term = "sx(mbmi)", main = "Mother body mass index", 
  xlab = "",ylab = "")
plot(zm1, term = "sx(mbmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6))
plot(zm1, term = "sx(mbmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6), rug = FALSE)
plot(zm1, term = "sx(mbmi)", main = "Mother body mass index", 
  xlab = "",ylab = "", ylim = c(-0.8, 0.6), rug = FALSE, 
  col.poly = NA, lwd = 1, lty = 1)
plot(zm1, term = "sx(mbmi)", main = "Mother body mass index", 
  xlab = "", ylab = "", ylim = c(-0.8, 0.6), rug = FALSE, 
  col.poly = NA, lwd = 1, lty = c(3, 1, 2, 2, 1))

plot(zm1, term = "sx(district):re", map = ZambiaBnd, pos = "topleft")
plot(zm1, term = "sx(district):re", map = ZambiaBnd, pos = "topleft")
plot(zm1, term = "sx(district):re", map = ZambiaBnd, pos = "topleft", 
  density = 30, swap = TRUE)
plot(zm1, term = "sx(district)", map = ZambiaBnd, pos = "topleft", 
  names = TRUE)
plot(zm1, term = "sx(district)", map = ZambiaBnd, pos = "topleft", 
  names = TRUE, cex.names = 0.8, cex.legend = 0.8)
plot(zm1, term = "sx(district)", map = ZambiaBnd, pos = "topleft", 
  range = c(-0.2, 0.2))
plot(zm1, term = "sx(district):mrf", map = ZambiaBnd, pos = "topleft", 
  c.select = "pcat95", at = c(-1, 0, 1), ncol = 3)
plot(zm1, term = "sx(district):mrf", map = ZambiaBnd, pos = "topleft", 
  c.select = "pcat80", at = c(-1, 0, 1), ncol = 3)

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
plot(zm1, term = "sx(district)", map = ZambiaBnd, pos = "topleft", 
  main = "Structured spatial effect", names = TRUE)
plot(zm1, term = "sx(district):re", map = ZambiaBnd, pos = "topleft", 
  density = 40, main = "Unstructured spatial effect", names = TRUE)
par(op)

## diagnostic plots
plot(zm1, which = 5:8, cex = 0.1)
plot(zm1, which = "intcpt-samples")

plot(zm1, term = c("memployment", "meducation"), 
  which = "coef-samples", ask = TRUE)

plot(zm1, term = "sx(mbmi)", which = "coef-samples")
plot(zm1, term = "sx(mbmi)", which = "coef-samples", acf = TRUE)

plot(zm1, term = "sx(mbmi)", which = "var-samples")
plot(zm1, term = "sx(mbmi)", which = "var-samples", acf = TRUE)
plot(zm1, term = "sx(mbmi)", which = "var-samples", acf = TRUE, lag.max = 400)

plot(zm1, term = "sx(agechild)", which = "var-samples")
plot(zm1, term = "sx(district)", which = "var-samples")
plot(zm1, term = "sx(district)", which = "var-samples")

## sensitivity analysis 
zm2 <- bayesx(stunting ~ memployment + meducation + urban + gender + 
  sx(mbmi, bs = "ps", a = 0.00001, b = 0.00001) +
  sx(agechild, bs = "ps", a = 0.00001, b = 0.00001) +
  sx(district, bs = "mrf", map = ZambiaBnd, a = 0.00001, b = 0.00001) +
  sx(district, bs = "re", a = 0.00001, b = 0.00001),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

zm3 <- bayesx(stunting ~ memployment + meducation + urban + gender +
  sx(mbmi, bs = "ps", a = 0.005, b = 0.005) +
  sx(agechild, bs = "ps", a = 0.005, b = 0.005) +
  sx(district, bs = "mrf", map = ZambiaBnd, a = 0.005, b = 0.005) +
  sx(district, bs = "re", a = 0.005, b = 0.005),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

zm4 <- bayesx(stunting ~ memployment + meducation + urban + gender +
  sx(mbmi, bs = "ps", a = 0.00005, b = 0.00005) +
  sx(agechild, bs = "ps", a = 0.00005, b = 0.00005) +
  sx(district, bs = "mrf", map = ZambiaBnd, a = 0.00005, b = 0.00005) +
  sx(district, bs = "re", a = 0.00005, b = 0.00005),
  iter = 12000, burnin = 2000, step = 10, data = ZambiaNutrition)

plot(c(zm1, zm2, zm3, zm4), term = "sx(mbmi)")
plot(c(zm1, zm2, zm3, zm4), term = "sx(agechild)")


## same model with REML
zm5 <- bayesx(stunting ~ memployment + meducation + urban + gender + 
  sx(mbmi) + sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) + sx(district, bs = "re"),
  method = "REML", data = ZambiaNutrition)
summary(zm5)
plot(zm5, map = ZambiaBnd)


## now use the stepwise algorithm
zm6 <- bayesx(stunting ~ memployment + meducation + urban + gender + 
  sx(mbmi) + sx(agechild) + sx(district, bs = "mrf", map = ZambiaBnd) + sx(district, bs = "re"),
  method = "STEP", data = ZambiaNutrition)
summary(zm6)
