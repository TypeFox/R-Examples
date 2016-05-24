options(warn = 2)

library(fit.models)
require(MASS)

################################################################################
# Test models with differenct formulas
################################################################################

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- c(gl(3,2), gl(3,1))
d.AD <- data.frame(treatment, outcome, counts)

trt <- glm(counts ~ treatment, family = poisson())
out <- glm(counts ~ outcome, family = poisson())

fm1 <- fit.models(trt, out)
print(fm1)
print(fm1.sum <- summary(fm1, correlation = TRUE))

pdf("fm1.pdf")
plot(fm1, 2)
plot(fm1, 3)
plot(fm1, 4)
plot(fm1, 5)
plot(fm1, 6)
plot(fm1, 7)
dev.off()

rm(fm1, fm1.sum)
unlink("fm1.pdf")


################################################################################
# Test models with different formulas
################################################################################


## an example with offsets from Venables & Ripley (2002, p.189)
utils::data(anorexia, package = "MASS")

offset <- glm(Postwt ~ Prewt + Treat + offset(Prewt), family = gaussian,
              data = anorexia)

no.offset <- glm(Postwt - Prewt ~ Treat, family = gaussian, data = anorexia)

fm2 <- fit.models(offset, no.offset)
print(fm2)
print(fm2.sum <- summary(fm2, correlation = TRUE))

pdf("fm2.pdf")
plot(fm2, 2)
plot(fm2, 3)
plot(fm2, 4)
plot(fm2, 5)
plot(fm2, 6)
plot(fm2, 7)
dev.off()

rm(fm2, fm2.sum)
unlink("fm2.pdf")


################################################################################
# Test models with different formulas and subsets
################################################################################

m1 <- glm(Postwt ~ Prewt + Treat, family = gaussian,
          data = subset(anorexia, Prewt > 75))

m2 <- glm(Postwt ~ Treat, family = gaussian, data = anorexia)

fm3 <- fit.models(m1, m2)
print(fm3)
print(fm3.sum <- summary(fm3, correlation = TRUE))

pdf("fm3.pdf")
plot(fm3, 2)
plot(fm3, 3)
plot(fm3, 4)
plot(fm3, 5)
plot(fm3, 6)
plot(fm3, 7)
dev.off()

rm(fm3, fm3.sum)
unlink("fm3.pdf")


################################################################################


