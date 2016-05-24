###################################################
### code chunk: Chap6init
###################################################
options(width=65, digits=5, show.signif.stars = FALSE)
date()
packageVersion("nlmeU")
packageVersion("nlme")
sessionInfo()

data(armd, package = "nlmeU")

###################################################
### code chunk: R6.1
###################################################
lm1.form <-                   # Fixed effects formula:(6.1)
    formula(visual ~ -1 + visual0 + time.f + treat.f:time.f )
vis.lm1.mf <- model.frame(lm1.form, armd)        # Model frame
vis.lm1.dm <- model.matrix(lm1.form, vis.lm1.mf) # Design matrix X
dim(vis.lm1.dm)               # Dimensions
(nms <- colnames(vis.lm1.dm)) # Long column names ...
nms <- abbreviate(nms)        # ... abbreviated
colnames(vis.lm1.dm) <- nms   # ... assigned.
head(vis.lm1.dm, n = 6)       # X matrix. Six rows.
attr(vis.lm1.dm, "contrasts") # Contrasts attribute.
contrasts(armd$treat.f)       # Contrasts for treat.f


###################################################
### code chunk: R6.2a
###################################################
lm6.1 <- lm(lm1.form, data = armd)         # M6.1:(6.1)
summ <- summary(lm6.1)                     # Summary
tT <- coef(summ)                           # beta, se, t-test
rownames(tT)                               # Fixed effects (beta) names
rownames(tT) <- abbreviate(rownames(tT))   # Abbreviated beta names
printCoefmat(tT, P.values = TRUE)
summ$sigma                                 # sigma


###################################################
### code chunk: R6.2b
###################################################
anova(lm6.1)                               # ANOVA table


###################################################
### code chunk: Simple syntax for Fig. 6.1a
###             Traditional graphics 
###################################################
plot(fitted(lm6.1), resid(lm6.1))          # Fig. 6.1a
abline(h = seq(-40, 40, by = 20), col = "grey")
abline(v = seq( 10, 80, by = 10), col = "grey")


###################################################
### code chunk: Elaborated syntax for Fig. 6.1a
###             Traditional graphics 
###################################################
ylim <- c(-50,50)
 plot(fitted(lm6.1), resid(lm6.1),
  ylim = ylim, type = "n", axes = FALSE, cex.lab = 0.9,
  ylab = "Residuals", 
  xlab = "Fitted values"
)  
abline(h = seq(-40,40, by = 20), col = "grey")
abline(v = seq(10, 80, by = 10), col = "grey")
points(fitted(lm6.1), resid(lm6.1))
axis(1, cex.axis = 0.9)
axis(2, cex.axis = 0.9)
box()

###################################################
### code chunk: Simple syntax for Fig. 6.1b
###             Traditional graphics 
###################################################
qqnorm(resid(lm6.1)); qqline(resid(lm6.1))  # Fig. 6.1b


###################################################
### code chunk: Elaborated syntax for Fig. 6.1b
###             Traditional graphics 
###################################################
qnDt <- qqnorm(resid(lm6.1),   
plot.it = FALSE
) 
plot(qnDt, 
  type = "n", 
  ylim = c(-50, 50), #  c(-3.8,3.8),  c(-50,50), 
  xlim = c(-3.8,3.8), #   c(-50,50), 
  main = "",                         # Normal Q-Q Plot
  ylab = "Sample Quantiles", 
  xlab = "Theoretical Quantiles", 
  axes = FALSE 
)
abline(h = seq(-40, 40, by = 20), col = "grey")
abline(v = -3:3, col = "grey")
qqline(resid(lm6.1))
points(qnDt)
axis(1, cex.axis = 0.9)
axis(2, cex.axis = 0.9)
box()


###################################################
### code chunk: R6.3
###################################################
require(nlme)              # Attach nlme package
fm6.1 <- gls(lm1.form,     # M6.1:(6.1)
             data = armd)
intervals(fm6.1)           # 95% CI for beta, sigma

###################################################
### code chunk: Syntax for Fig. 6.1
###################################################
plot(predict(fm6.1), residuals(fm6.1))   # Same as Fig. 6.1a
qqnorm(residuals(fm6.1))                 # Same as Fig. 6.1b
qqline(residuals(fm6.1))

### sessionInfo()
sessionInfo()
detach(package:nlme)

