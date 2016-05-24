qbeta(c(0.025,0.975),38+1,62+1)
binom.test(38,100)$conf.int            # for comparison
prop.test(38,100)$conf.int             # for comparison
