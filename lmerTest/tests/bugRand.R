require(lmerTest)
load(system.file("testdata","cltlike.RData", package="lmerTest"))

cltlike$datef <- as.factor(cltlike$date)

fm2b <- lmer(liking ~ pig.type + pig.type:landro + pig.type:lskatol +
                             sequence +
                             (1|pig) +
                             (1|session) +
                             (1|serving) +
                             (1|datef) +
                             (1|consumer),
                           data=cltlike)

VarCorr(fm2b)

rnd <- rand(fm2b)
TOL <- 1e-4
stopifnot(all.equal(rnd$rand.table[,"p.value"], c(0.02952, 0.7847, 1,
                                                1, 1.509e-14), tol = TOL), TRUE)


m.carrots <- lmer(Preference ~ sens1 + sens2 +
                    (1 + sens1 + sens2|Consumer) + (1|product), data = carrots)

rnd <- rand(m.carrots)

TOL <- 1e-4
stopifnot(all.equal(rnd$rand.table[,"p.value"], c(0.609,0.02284, 2.318e-05), tol = TOL), TRUE)


## check that step and rand functions give same results


m.carrots2 <- lmer(Preference ~ sens2*Homesize*Age +
                    (1 + sens2|Consumer) + (1|product), data = carrots)

rnd.car <- rand(m.carrots2)



st.car <- step(m.carrots2)


stopifnot(all.equal(rnd.car$rand.table$p.value, st.car$rand.table$p.value))


