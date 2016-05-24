require(Cprob)

op <- options()
options(warn = 1)

### test 1
data(mgus, package = "Cprob")

aa <- cpf(Hist(time, ev) ~ 1, mgus)

aa

summary(aa)

### test 2
mgus$A <- ifelse(mgus$age < 64, 0, 1)

bb <- cpf(Hist(time, ev) ~ A, mgus)

bb

summary(bb)

all.equal(bb[]$cp, bb$cp)
all.equal(bb[c(1, 2)]$cp, bb$cp)

all.equal(bb[1]$cp, bb$cp[1:bb$size.strata[1]])

### test 3
fit <- cpfpo(Hist(time, ev) ~ age + creat, mgus,
             tis=seq(10, 30, 0.3), w=rep(1,67))

### test 4
time <- c(rep(2, 20), rep(3, 60), rep(4, 20))
ev <- c(rep(1, 20), rep(2, 60), rep(0, 20))
data <- data.frame(time, ev)

a <- predict(cpf(Hist(time, ev) ~ 1, data), 3)$cp
all.equal(a, 1/2)

### test 5
cutoffs <- quantile(mgus$time, probs = seq(0, 1, 0.05))[-1]

fit1 <- pseudocpf(Hist(time, ev) ~ age + creat, mgus, id = id, timep = cutoffs,
                 corstr = "independence", scale.value = TRUE)

fit2 <- pseudocpf(Hist(time, ev) ~ age + creat, mgus, id = id, timep = cutoffs,
                 corstr = "independence", scale.value = TRUE, jack = TRUE)

fit1
fit2

summary(fit1)
summary(fit2)


### test 6: playing with cens.code
mm <- mgus
mm$ev <- ifelse(mm$ev == 0, 4, mm$ev)

test1 <- pseudocpf(Hist(time, ev, cens.code = "4") ~ age + creat, mm, id = id, timep = cutoffs,
                   corstr = "independence", scale.value = TRUE)

test2 <- cpfpo(Hist(time, ev, cens.code = "4") ~ age + creat, mm,
               tis=seq(10, 30, 0.3), w=rep(1,67))

all(test2$alpha == fit$alpha)

all(summary(test1)$coef == summary(fit1)$coef)

options(op)
