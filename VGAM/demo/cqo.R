# Demo for constrained quadratic ordination (CQO; aka
# canonical Gaussian ordination)


data(hspider, package = "VGAM")
hspider[, 1:6] <- scale(hspider[, 1:6])  # standardize environmental vars


## Rank-1 model (unequal tolerances, deviance = 1176.0)

set.seed(123)
p1 <-
  cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
            Pardlugu, Pardmont, Pardnigr, Pardpull, Trocterr, Zoraspin) ~
      WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux, 
      quasipoissonff, data = hspider, 
      Bestof = 10, Crow1positive = FALSE, eq.tolerances = FALSE,
      I.tolerances = FALSE)

par(mfrow = c(3, 3))

lvplot(p1, lcol = 1:12, llwd = 2, llty = 1:12, y = TRUE, pch = 1:12,
       pcol = 1:12, las = 1, main = "Hunting spider data")

print(cancoef(p1), digits = 3)
print(Coef(p1), digits = 3)

# trajectory plot
trplot(p1, which = 1:3, log = "xy", type = "b", lty = 1, 
       col = c("blue", "orange", "green"), lwd = 2, label = TRUE) -> ii
legend(0.00005, 0.3, paste(ii$species[, 1], ii$species[, 2], sep = " and "),
       lwd = 2, lty = 1, col = c("blue", "orange", "green"))
abline(a = 0, b = 1, lty = "dashed")




## Rank-2 model (equal tolerances, deviance = 856.5) 

set.seed(111)
r2 <-
  cqo(cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi, Auloalbi,
            Pardmont, Pardnigr, Pardpull, Trocterr) ~
      WaterCon + BareSand + FallTwig + CoveMoss + CoveHerb + ReflLux,
      quasipoissonff, data = hspider, Rank = 2,
      Bestof = 10, I.tolerances = TRUE,
      eq.tolerances = TRUE, Crow1positive = c(FALSE, FALSE))
print(ccoef(r2), digits = 3)
print(Coef(r2), digits = 3)

clr <- (1:(10+1))[-7]  # Omit yellow colour
adj <- c(-0.1, -0.1, -0.1, 1.1, 1.1, 1.1, -0.1, -0.1, -0.1, 1.1)
# With C arrows
lvplot(r2, label = TRUE, xlim = c(-2.8, 5.0), ellipse = FALSE, C = TRUE,
       Cadj = c(1.1, -0.1, 1.2, 1.1, 1.1, -0.1), adj = adj,
       las = 1, chull = TRUE, pch = "+", pcol = clr, sites = TRUE)

# With circular contours
lvplot(r2, label = TRUE, xlim = c(-2.8, 5.0), ellipse = TRUE, C = FALSE,
       Cadj = c(1.1, -0.1, 1.2, 1.1, 1.1, -0.1), adj = adj,
       las = 1, chull = TRUE, pch = "+", pcol = clr, sites = TRUE)

# With neither C arrows or circular contours
lvplot(r2, label = TRUE, xlim = c(-2.8, 5.0), ellipse = FALSE, C = FALSE,
       Cadj = c(1.1, -0.1, 1.2, 1.1, 1.1, -0.1), adj = adj,
       las = 1, chull = TRUE, pch = "+", pcol = clr, sites = TRUE)

# Perspective plot 
persp(r2, xlim = c(-5, 5), ylim = c(-3, 6), theta = 50, phi = 20)



## Gaussian logit regression
## Not recommended actually because the number of sites is far too low.
## Deviance = 154.6, equal tolerances.

ybin <- with(hspider,
             0 + (cbind(Alopacce, Alopcune, Alopfabr, Arctlute, Arctperi,
                        Auloalbi, Pardlugu, Pardmont, Pardnigr, Pardpull,
                        Trocterr, Zoraspin) > 0))  # Matrix of 0s and 1s
colnames(ybin) <- paste0(colnames(ybin), ".01")
hspider <- data.frame(hspider, ybin)

set.seed(1312)
b1 <- cqo(ybin[, -c(1, 5)] ~ WaterCon + BareSand + FallTwig + CoveMoss +
         CoveHerb + ReflLux, quasibinomialff(mv = TRUE), 
         Bestof = 4, I.tolerances = TRUE, 
         data = hspider, eq.tolerances = TRUE, Crow1positive = FALSE)
lvplot(b1, type = "predictors", llwd = 2, las = 1, ylab = "logit mu",
       ylim = c(-20, 11), lcol = 1:10)
c1 <- Coef(b1)
cts <- c("Trocterr", "Pardmont", "Alopfabr", "Arctlute")
text(c1@Optimum[1, cts], logit(c1@Maximum[cts])+1.0, cts)

round(t(Coef(b1, I.tolerances = FALSE)@C), dig = 3)

# On the probability scale
lvplot(b1, type = "fitted", llwd = 2, las = 1, llty = 1,
       ylab = "Probability of presence",
       ylim = c(0, 1), lcol = 1:10)


