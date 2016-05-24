
if (requireNamespace("car") && require("effects")){
  
  data(Duncan, package="car")

  mi <- with(Duncan, mean(income))
  me <- with(Duncan, mean(education))
  med <- with(Duncan, median(education))
  
  # (1) focal: factor, constant: polynomial
  
  mod.1 <- lm(prestige ~ type + poly(income, degree=2, raw=TRUE), data=Duncan)
  
  X <- matrix(c(1, 0, 0, mi, mi^2,
                1, 1,  0, mi, mi^2,
                1, 0, 1, mi, mi^2),
              nrow=3, ncol=5, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(matrix(X %*% coef(mod.1))), as.vector(Effect("type", mod.1)$fit)))) 
    stop("failed Test 1-1")
 
  
  # (2) focal: polynomial, constant: factor
  
  X <- matrix(c(1, 0.4, 2/15, 10, 10^2,
                  1, 0.4, 2/15, 40, 40^2,
                  1, 0.4, 2/15, 70, 70^2),
                nrow=3, ncol=5, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(Effect("income", mod.1, xlevels=list(income=c(10, 40, 70)))$fit), 
            as.vector(matrix(X %*% coef(mod.1))))))
    stop("failed test 1-2")
  
  # (2a) As in (2), but without specifying xlevels
  X <- matrix(c(1, 0.4, 2/15, 10, 10^2,
                1, 0.4, 2/15, 20, 20^2,
                1, 0.4, 2/15, 30, 30^2,
                1, 0.4, 2/15, 40, 40^2,
                1, 0.4, 2/15, 50, 50^2,
                1, 0.4, 2/15, 60, 60^2,
                1, 0.4, 2/15, 70, 70^2,
                1, 0.4, 2/15, 80, 80^2),
              nrow=8, ncol=5, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(Effect("income", mod.1)$fit), 
                        as.vector(matrix(X %*% coef(mod.1)))))) 
    stop("failed test 1-2a")
  
  # (3) focal: factor*polynomial, constant: polynomial
  
  mod.2 <- lm(prestige ~ type*poly(income, degree=2, raw=TRUE) + 
                poly(education, degree=2, raw=TRUE), data=Duncan)
  
  X <- matrix(c(1, 0, 0, 10, 10^2, me, me^2, 0, 0, 0, 0,
                1, 1, 0, 10, 10^2, me, me^2, 10, 0, 10^2, 0,
                1, 0, 1, 10, 10^2, me, me^2, 0, 10, 0,  10^2,
                1, 0, 0, 70, 70^2, me, me^2, 0, 0, 0, 0,
                1, 1, 0, 70, 70^2, me, me^2, 70, 0, 70^2, 0,
                1, 0, 1, 70, 70^2, me, me^2, 0, 70, 0,  70^2),
              nrow=6, ncol=11, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(Effect(c("type", "income"), mod.2, xlevels=list(income=c(10, 70)))$fit),
            as.vector(matrix(X %*% coef(mod.2), 3, 2)))))
    stop("failed test 1-3")
  
  # (4) focal: polynomial, constant: factor*polynomial
  
  X <- matrix(c(1, 0.4, 2/15, mi, mi^2, 10, 10^2, 0.4*mi, 2/15*mi, 0.4*mi^2, 2/15*mi^2,
                1, 0.4, 2/15, mi, mi^2, 40, 40^2, 0.4*mi, 2/15*mi, 0.4*mi^2, 2/15*mi^2,
                1, 0.4, 2/15, mi, mi^2, 70, 70^2, 0.4*mi, 2/15*mi, 0.4*mi^2, 2/15*mi^2),
              nrow=3, ncol=11, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(Effect("education", mod.2, xlevels=list(education=c(10, 40, 70)))$fit),
            as.vector(X %*% coef(mod.2)))))
    stop("failed test 1-4")
  
  # (5) repeat of (3) with medians rather than means
  
  X <- matrix(c(1, 0, 0, 10, 10^2, med, med^2, 0, 0, 0, 0,
                1, 1, 0, 10, 10^2, med, med^2, 10, 0, 10^2, 0,
                1, 0, 1, 10, 10^2, med, med^2, 0, 10, 0,  10^2,
                1, 0, 0, 70, 70^2, med, med^2, 0, 0, 0, 0,
                1, 1, 0, 70, 70^2, med, med^2, 70, 0, 70^2, 0,
                1, 0, 1, 70, 70^2, med, med^2, 0, 70, 0,  70^2),
              nrow=6, ncol=11, byrow=TRUE)
  if (!isTRUE(all.equal(as.vector(Effect(c("type", "income"), mod.2, 
                             xlevels=list(income=c(10, 70)), typical=median)$fit),
            as.vector(X %*% coef(mod.2)))))
    stop("failed test 1-5")
  
  # (6) focal: factor*polynomial, constant: polynomial, using predict() & orthog. polys.
  
  mod.3 <- lm(prestige ~ type*poly(income, degree=2) + poly(education, degree=2), data=Duncan)
  
  if (!isTRUE(all.equal(as.vector(predict(mod.3, 
                              newdata=data.frame(income=c(10, 10, 10, 70, 70, 70), 
                                                type=factor(c("bc", "prof", "wc", "bc", "prof", "wc")),
                                                education=mean(Duncan$education)))),
            as.vector(Effect(c("type", "income"), mod.3, xlevels=list(income=c(10, 70)))$fit))))
    stop("failed test 1-6")
  
  # (7) focal: factor, constant: poly*poly
  
  mod.4 <- lm(prestige ~ type + poly(income, 2)*poly(education, 2), data=Duncan)
  
  if (!isTRUE(all.equal(as.vector(Effect("type", mod.4)$fit),
            as.vector(predict(mod.4, newdata=data.frame(type=c("bc", "prof", "wc"), 
                                                        income=rep(mi, 3), education=rep(me, 3)))))))
    stop("failed test 1-7")
  
  # (8) focal: factor, constant: 2nd deg polynomial in 2 Xs
  
  mod.5 <- lm(prestige ~ type + poly(income, education, degree=2), data=Duncan)
  
  if (!isTRUE(all.equal(as.vector(Effect("type", mod.5)$fit),
            as.vector(predict(mod.5, newdata=data.frame(type=c("bc", "prof", "wc"), 
                                    income=rep(mi, 3), education=rep(me, 3)))))))
    stop("failed test 1-8")
  
  # (9) focal: covariate, constant: 2 factors and 1 covariate, 3-way interaction
  
  data(Mroz, package="car")
  mod.6 <- lm(lwg ~ inc + age*hc*wc, data=Mroz)
  mage <- with(Mroz, mean(age))
  mhc <- with(Mroz, mean(hc == "yes"))
  mwc <- with(Mroz, mean(wc == "yes"))
  hc <- rep(mhc, 3)
  wc <- rep(mwc, 3)
  age <- rep(mage, 3)
  X <- cbind(1, c(10, 40, 80), age, hc, wc, age*hc, age*wc, hc*wc, age*hc*wc)
  if (!isTRUE(all.equal(as.vector(Effect("inc", mod.6, xlevels=list(inc=c(10, 40, 80)))$fit),
                        as.vector(X %*% coef(mod.6)))))
      stop("failed test 1-8")
}
