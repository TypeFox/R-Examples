options(warn = 2)

data(stackloss)

library(fit.models)
require(MASS)

# add rlm to the lmfm fit.models class
fmclass.add.class("lmfm", "rlm")

################################################################################
# Test original syntax
################################################################################

fm1 <- fit.models(list(Robust = "rlm", LS = "lm"), stack.loss ~ ., data = stackloss)
print(fm1)
print(fm1.sum <- summary(fm1, correlation = TRUE))

pdf("fm1.pdf")
plot(fm1, 2)
plot(fm1, 3)
plot(fm1, 4)
plot(fm1, 5)
plot(fm1, 6)
plot(fm1, 7)
plot(fm1, 8)
plot(fm1, 9)
plot(fm1, 10)
dev.off()

rm(fm1, fm1.sum)
unlink("fm1.pdf")


################################################################################
# Test models with different subsets
################################################################################

complete <- lm(stack.loss ~ ., data = stackloss)
clean <- lm(stack.loss ~ ., data = stackloss, subset = -c(1, 2, 4, 21))
fm2 <- fit.models(Clean = clean, Complete = complete)
print(fm2)
print(fm2.sum <- summary(fm2, correlation = TRUE))

pdf("fm2.pdf")
plot(fm2, 2)
plot(fm2, 3)
plot(fm2, 4)
plot(fm2, 5)
plot(fm2, 6)
plot(fm2, 7)
plot(fm2, 8)
plot(fm2, 9)
plot(fm2, 10)
dev.off()

rm(fm2, fm2.sum)
unlink("fm2.pdf")


################################################################################
# Test models with different formulas
################################################################################

m1 <- lm(stack.loss ~ Air.Flow + Water.Temp, data = stackloss)
m2 <- lm(stack.loss ~ Water.Temp + Acid.Conc., data = stackloss)
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
plot(fm3, 8)
plot(fm3, 9)
plot(fm3, 10)
dev.off()

rm(fm3, fm3.sum)
unlink("fm3.pdf")


################################################################################
# Test models with different formulas and different subsets
################################################################################

m1 <- lm(stack.loss ~ Air.Flow + Water.Temp, data = stackloss, subset = -c(1, 2, 4, 21))
m2 <- lm(stack.loss ~ Water.Temp + Acid.Conc., data = stackloss)
fm4 <- fit.models(m1, m2)
print(fm4)
print(fm4.sum <- summary(fm4, correlation = TRUE))

pdf("fm4.pdf")
plot(fm4, 2)
plot(fm4, 3)
plot(fm4, 4)
plot(fm4, 5)
plot(fm4, 6)
plot(fm4, 7)
plot(fm4, 8)
plot(fm4, 9)
plot(fm4, 10)
dev.off()

rm(fm4, fm4.sum)
unlink("fm4.pdf")


################################################################################
# Test simple linear regression
################################################################################

fm5 <- fit.models(list(Robust = "rlm", LS = "lm"), stack.loss ~ Acid.Conc., data = stackloss)
print(fm5)
print(fm5.sum <- summary(fm5, correlation = TRUE))

pdf("fm5.pdf")
plot(fm5, 2)
plot(fm5, 3)
plot(fm5, 4)
plot(fm5, 5)
plot(fm5, 6)
plot(fm5, 7)
plot(fm5, 8)
plot(fm5, 9)
plot(fm5, 10)
plot(fm5, 11)
dev.off()

rm(fm5, fm5.sum)
unlink("fm5.pdf")


################################################################################


