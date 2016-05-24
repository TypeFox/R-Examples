require(lmerTest)

### example from the paper GLM SAS 101 report
a <- factor(c(1,1,1,2,2,2,2,2,1,2))
b <- factor(c(1,1,2,1,2,2,2,2,2,1))
f=factor(c(1,2,1,2,1,2,1,2,1,2))
y <- c(12,14,11,20,17,23,35,46,15,16)
dd <- data.frame(a=a, b=b, y=y, f=f)

## check type 2 is order independent
model <- lmer(y ~ a*b + (1|f), data=dd)
model2 <- lmer(y ~ b*a + (1|f), data=dd)
an <- anova(model, type=2)
an2 <- anova(model2, type=2)
stopifnot(all.equal(an,an2[c(2,1,3),], check.names = FALSE, check.attributes = FALSE))

## check the results are the same as from SAS proc mixed
stopifnot(all.equal(an[,"F.value"], c(3.90131, 1.32753, 0.99565), tol = 1e-5))


m.carrots <- lmer(Preference ~ sens2*Homesize
                  +(1+sens2|Consumer), data=carrots)
ancar <- anova(m.carrots, type=2)

stopifnot(all.equal(ancar[,"F.value"], c(54.8361, 5.16138, 1.03035), tol = 1e-3))

m <- lmer(Informed.liking ~ Product*Age 
          + (1|Consumer) , data=ham)
an <- anova(m, type=2)

stopifnot(all.equal(an[,"F.value"], c(2.48135, .005387, 1.48451), tol = 1e-5))



fm <- lmer(Preference ~ sens2*Homesize*sens1 + (1|product),
              data=carrots)

#check type2
ant2 <- anova(fm, type=2)

stopifnot(all.equal(ant2[,"F.value"], c(16.4842, 14.0010, .526076, 1.18144,
                                        .107570, .335177, 1.05946), tol = 1e-4))

ant3 <- anova(fm, type=3)

stopifnot(all.equal(ant3[,"F.value"], c(16.9140, 14.0010,.481148, 1.18144,
                                        .074201, .335177, 1.05946), tol = 1e-4))
