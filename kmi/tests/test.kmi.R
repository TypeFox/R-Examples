require(kmi)

### test 1 

id <- 1:10
time <- 1:10
ev <- c(0, 0, 0, rep(1, 3), rep(2, 4))
data <- data.frame(id, time, ev)

aa <- kmi(Surv(time, ev != 0) ~ 1, data, etype = ev)

bb <- aa$imputed.data
test1 <- lapply(bb, function(x) {
    x[aa$original.data$id %in% 7:10, ]
})
test1 <- do.call(rbind, test1)

all(test1[, 1] == 11)
all(test1[, 2] == 0)

aa <- kmi(Surv(time, ev != 0) ~ 1, data, etype = ev, failcode = 2)

bb <- aa$imputed.data
test1 <- lapply(bb, function(x) {
    x[aa$original.data$id %in% 4:6, ]
})
test1 <- do.call(rbind, test1)

all(test1[, 1] == 11)
all(test1[, 2] == 0)

aa <- kmi(Surv(time, ev != 0) ~ 1, data, etype = ev, epsilon = 2)

bb <- aa$imputed.data
test1 <- lapply(bb, function(x) {
    x[aa$original.data$id %in% 7:10, ]
})
test1 <- do.call(rbind, test1)

all(test1[, 1] == 12)
all(test1[, 2] == 0)

aa <- kmi(Surv(time, ev != 0) ~ 1, data, etype = ev, nimp = 13)
length(aa$imputed.data) == 13


### test 2

set.seed(198740)
time <- rexp(100)
ev <- sample(c(0, 1, 2), 100, replace = TRUE)
cov <- rbinom(100, 1, 0.5)
dd <- data.frame(time, ev, cov)

## add a test when etype is a factor
dd$status <- factor(ifelse(dd$ev == 0, "cens", ifelse(dd$ev == 1, "rel", "dc")))

set.seed(1440293)
dat.kmi <- kmi(Surv(time, ev != 0) ~ 1, dd, etype = ev, nimp = 5)
set.seed(1440293)
dat.kmi.fact <- kmi(Surv(time, status != "cens") ~ 1, dd, etype = status,
                    nimp = 5, failcode = "rel")
set.seed(1440293)
dat.kmi.mixed1 <- kmi(Surv(time, ev != 0) ~ 1, dd, etype = status,
                      nimp = 5, failcode = "rel")
set.seed(1440293)
dat.kmi.mixed2 <- kmi(Surv(time, status != "cens") ~ 1, dd, etype = ev,
                      nimp = 5, failcode = 1)


fit.kmi <- cox.kmi(Surv(time, ev == 1) ~ cov, dat.kmi)
fit.kmi.fact <- cox.kmi(Surv(time, status == "rel") ~ cov, dat.kmi)
fit.kmi.mixed1 <- cox.kmi(Surv(time, status == "rel") ~ cov, dat.kmi.mixed1)
fit.kmi.mixed2 <- cox.kmi(Surv(time, ev == 1) ~ cov, dat.kmi.mixed2)

all.equal(coef(fit.kmi), coef(fit.kmi.fact))
all.equal(coef(fit.kmi), coef(fit.kmi.mixed1))
all.equal(coef(fit.kmi), coef(fit.kmi.mixed2))
all.equal(fit.kmi$variance, fit.kmi.fact$variance)
all.equal(fit.kmi$variance, fit.kmi.mixed1$variance)
all.equal(fit.kmi$variance, fit.kmi.mixed2$variance)

## fit.kmi

## summary(fit.kmi)

## fit.kmi.fact

## summary(fit.kmi.fact)

## avec bootstrap

set.seed(867988)
dat.kmib <- kmi(Surv(time, ev != 0) ~ 1, dd, etype = ev, nimp = 5,
                boot = TRUE, nboot = 5)
set.seed(867988)
dat.kmib.fact <- kmi(Surv(time, ev != 0) ~ 1, dd, etype = status, nimp = 5,
                     boot = TRUE, nboot = 5, failcode = "rel")

fit.kmib <- cox.kmi(Surv(time, ev == 1) ~ cov, dat.kmib)
fit.kmib.fact <- cox.kmi(Surv(time, status == "rel") ~ cov, dat.kmib.fact)

all.equal(coef(fit.kmib), coef(fit.kmib.fact))

fit.kmib
fit.kmib.fact

summary(fit.kmib)


### test 4

data(icu.pneu)

set.seed(1313)
dat <- kmi(Surv(start, stop, status) ~ 1, data = icu.pneu,
           etype = event, id= id, failcode = 2, nimp = 5)

### add a factor for testing purposes
icu.pneu$ev.fact <- factor(ifelse(icu.pneu$event == 3, "disch", "death"))

icu.pneu$ev <- icu.pneu$event
icu.pneu$ev[icu.pneu$status == 0] <- 0

set.seed(1313)
dat2 <- kmi(Surv(start, stop, ev != 0) ~ 1, data = icu.pneu,
           etype = ev, id= id, failcode = 2, nimp = 5)

set.seed(1313)
dat3 <- kmi(Surv(start, stop, status) ~ 1, data = icu.pneu,
            etype = ev.fact, id = id, failcode = "death", nimp = 5)

a <- logical(5)
for (i in 1:5) a[i] <- all.equal(dat$imputed.data[[i]][, 1], dat2$imputed.data[[i]][, 1])
a

fit.kmi <- cox.kmi(Surv(start, stop, event == 2) ~ pneu, dat)

fit.kmi2 <- cox.kmi(Surv(start, stop, ev == 2) ~ pneu, dat2)

fit.kmi3 <- cox.kmi(Surv(start, stop, ev.fact == "death") ~ pneu, dat3)

all.equal(fit.kmi$coefficients, fit.kmi2$coefficients)
all.equal(coef(fit.kmi), coef(fit.kmi3))
all.equal(fit.kmi$variance, fit.kmi2$variance)

fit.kmi

fit.kmi2

fit.kmi3

## avec bootstrap

set.seed(598085)
dat <- kmi(Surv(start, stop, status) ~ 1, data = icu.pneu,
           etype = event, id= id, failcode = 2, nimp = 5,
           boot = TRUE, nboot = 5)

set.seed(598085)
dat2 <- kmi(Surv(start, stop, ev != 0) ~ 1, data = icu.pneu,
            etype = ev, id= id, failcode = 2, nimp = 5,
            boot = TRUE, nboot = 5)

set.seed(598085)
dat3 <- kmi(Surv(start, stop, status) ~ 1, data = icu.pneu,
            etype = ev.fact, id = id, failcode = "death", nimp = 5,
            boot = TRUE, nboot = 5)


a <- logical(5)
for (i in 1:5) a[i] <- all.equal(dat$imputed.data[[i]][, 1], dat2$imputed.data[[i]][, 1])
a

fit.kmi <- cox.kmi(Surv(start, stop, event == 2) ~ pneu, dat)

fit.kmi2 <- cox.kmi(Surv(start, stop, ev == 2) ~ pneu, dat2)

fit.kmi3 <- cox.kmi(Surv(start, stop, ev.fact == "death") ~ pneu, dat3)

all.equal(fit.kmi$coefficients, fit.kmi2$coefficients)
all.equal(coef(fit.kmi), coef(fit.kmi3))
all.equal(fit.kmi$variance, fit.kmi2$variance)

fit.kmi

fit.kmi2

fit.kmi3

### with covariates
## classic
set.seed(1)
dd$juhu <- rnorm(nrow(dd))

set.seed(78223)
imp.dd <- kmi(Surv(time, ev != 0) ~ juhu, dd,
              etype = ev, nimp = 5)
set.seed(44889)
imp.ddb <- kmi(Surv(time, ev != 0) ~ juhu, dd, nimp = 5,
              etype = ev, boot = TRUE, nboot = 5)

summary(cox.kmi(Surv(time, ev == 1) ~ cov, imp.dd))
summary(cox.kmi(Surv(time, ev == 1) ~ cov, imp.ddb))

## time-dependent covariates
set.seed(9763)
imp.dat.c <- kmi(Surv(start, stop, status) ~ age + sex,
                 data = icu.pneu, etype = event, id = id,
                 failcode = 2, nimp = 5)

set.seed(19832)
imp.dat.cb <- kmi(Surv(start, stop, status) ~ age + sex,
                  data = icu.pneu, etype = event, id = id,
                  failcode = 2, nimp = 5, boot = TRUE,
                  nboot = 5)

summary(cox.kmi(Surv(start, stop, event == 2) ~ pneu,
                imp.dat.c))
summary(cox.kmi(Surv(start, stop, event == 2) ~ pneu,
                imp.dat.cb))
