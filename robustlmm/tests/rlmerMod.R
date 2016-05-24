## Test disabled
quit()

require(robustlmm)

b <- robustlmm:::b
b.rlmerMod <- robustlmm:::b.rlmerMod
u <- robustlmm:::u
u.rlmerMod <- robustlmm:::u.rlmerMod
std.b <- robustlmm:::std.b
Lambda <- robustlmm:::Lambda

## one-way anova
fm1 <- lmer(Yield ~ (1 | Batch), Dyestuff)
rfm1 <- rlmer(Yield ~ (1 | Batch), Dyestuff, doFit = FALSE)

## correlated random effects
fm4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
rfm4 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit = FALSE)

sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])

## correlated random effects and zeroes in the variance components
fm5 <- lmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)
rfm5 <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2, doFit = FALSE)


## test updating fixef
robustlmm:::fixef(rfm1) <- 1500
stopifnot(fixef(rfm1) == 1500,
          all(rfm1@resp$mu == 1500 + robustlmm:::getZ(rfm1) %*% b(rfm1)),
          all(rfm1@resp$wtres == rfm1@resp$y - rfm1@resp$mu)## ,
          ## all(resid(rfm1) == rfm1@resp$res * rfm1@resp$weights)
          )
robustlmm:::fixef(rfm4) <- c(250, 10)
stopifnot(fixef(rfm4) == c(250, 10),
          all(rfm4@resp$mu == robustlmm:::getX(rfm4) %*% c(250, 10) +
              robustlmm:::getZ(rfm4) %*% b(rfm4)),
          all(rfm4@resp$wtres == rfm4@resp$y - rfm4@resp$mu)## ,
          ## all(resid(rfm4) == rfm4@resp$res * rfm4@resp$weights)
          )

## test updating theta
robustlmm:::theta(rfm1, fit.effects = TRUE, update.sigma = FALSE) <- 2
stopifnot(theta(rfm1) == 2,
          diag(Lambda(rfm1)) == 2)
robustlmm:::theta(rfm4, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 0, 2)
stopifnot(theta(rfm4) == c(1, 0, 2),
          diag(Lambda(rfm4)) == c(1, 2))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 0, 0, 0)
stopifnot(theta(rfm5) == c(1, 0, 0, 0))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 1, 0, 0)
stopifnot(theta(rfm5) == c(1, 1, 0, 0))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(0, 1, 0, 1)
stopifnot(theta(rfm5) == c(0, 0, 0, 1))

## test updating u
## test integrity of rfm1:
stopifnot(all(Lambda(rfm1) %*% u(rfm1) == b(rfm1)),
          all(std.b(rfm1, 1, Matrix(b(rfm1))) == u(rfm1)))
rfm1a <- rfm1
rfm1b <- rfm1
robustlmm:::theta(rfm1a, fit.effects = TRUE, update.sigma = FALSE) <- 3
robustlmm:::theta(rfm1b, fit.effects = TRUE, update.sigma = FALSE) <- 3
## before setting u or b
stopifnot(all.equal(u(rfm1a), u(rfm1b)),
          all.equal(drop(Lambda(rfm1a) %*% u(rfm1a)), unname(b(rfm1a))),
          all.equal(std.b(rfm1a, 1, Matrix(b(rfm1a))), unname(u(rfm1a))),
          all.equal(b(rfm1a), b(rfm1b)),
          all.equal(drop(Lambda(rfm1b) %*% u(rfm1b)), unname(b(rfm1b))),
          all.equal(std.b(rfm1b, 1, Matrix(b(rfm1b))), unname(u(rfm1b))))
robustlmm:::u(rfm1a) <- u(rfm1)
robustlmm:::b(rfm1b) <- drop(Lambda(rfm1b) %*% u(rfm1))
## after:
stopifnot(all.equal(u(rfm1a), u(rfm1)),
          all.equal(unname(b(rfm1b)), drop(Lambda(rfm1b) %*% u(rfm1))),
          all.equal(u(rfm1a), u(rfm1b)),
          all.equal(drop(Lambda(rfm1a) %*% u(rfm1a)), unname(b(rfm1a))),
          all.equal(std.b(rfm1a, 1, Matrix(b(rfm1a))), unname(u(rfm1a))),
          all.equal(b(rfm1a), b(rfm1b)),
          all.equal(drop(Lambda(rfm1b) %*% u(rfm1b)), unname(b(rfm1b))),
          all.equal(std.b(rfm1b, 1, Matrix(b(rfm1b))), unname(u(rfm1b))))

## test dependency on initial values of theta
testInit <- function(formula, data, ...) {
    fm_1     <- lmerNoFit(formula, data)
    fm_10000 <- lmerNoFit(formula, data, initTheta = 10000)
    fm_0     <- lmerNoFit(formula, data, initTheta = 0)
    
    o1 <- capture.output(print(rlmer(formula, data, ..., init = fm_1)))
    o2 <- capture.output(print(rlmer(formula, data, ..., init = fm_10000)))
    o3 <- capture.output(print(rlmer(formula, data, ..., init = fm_0)))
    
    stopifnot(all.equal(o1, o2),
              all.equal(o1, o3))
}

testInit(Yield ~ (1 | Batch), Dyestuff)
testInit(Yield ~ (1 | Batch), Dyestuff, rho.e = smoothPsi, rho.b = smoothPsi)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
