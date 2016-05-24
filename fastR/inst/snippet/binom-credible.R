qbeta(c(0.025,0.975),20+1,30+1)
binom.test(20,50)$conf.int             # for comparison
prop.test(20,50)$conf.int              # for comparison
