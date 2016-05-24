#
# Tests of integer case weights
#
library(deming)

data1 <- ferritin[ferritin$period==2,]
data1$wt <- c(1:10, 10:1)
data2 <- data1[rep(1:20, c(1:10, 10:1)),]

fit1a <- deming(new.lot ~ old.lot, data1, weight=wt)
fit1b <- deming(new.lot ~ old.lot, data2)
all.equal(coef(fit1a), coef(fit1b))

fit2a <- deming(new.lot ~ old.lot, data1, cv=TRUE, weight=wt, conf=0)
fit2b <- deming(new.lot ~ old.lot, data2, cv=TRUE, conf=0)
all.equal(coef(fit2a), coef(fit2b))

# Variable case wieghts do not replicate exactly for Thiel-Sen regression
#  It is a problem of discreteness. 
# However making all case weights the same gives a valid test
data3 <- data1[rep(1:20, 2),]
data1$wt2 <- rep(2,20)
for (sym in c(FALSE, TRUE)) {
    fit3a <- thielsen(new.lot ~ old.lot, data1, weight=wt2, conf=0, 
                      symmetric=sym)
    fit3b <- thielsen(new.lot ~ old.lot, data3, conf=0, symmetric=sym)
    print(all.equal(coef(fit3a), coef(fit3b)))
}
 
for (m in 1:3) {
     fit4a <- pbreg(new.lot ~ old.lot, data1, weight=wt2, conf=0,
                    method=m)
     fit4b <- pbreg(new.lot ~ old.lot, data3, conf=0, method=m)
     print(all.equal(coef(fit4a), coef(fit4b)))
 }

