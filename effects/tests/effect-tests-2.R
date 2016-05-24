
if (requireNamespace("car") && require("effects")){

  # plots should show fitted values directly on plotted effect, and must be checked visually
  # numbering corresponds to effect-test-1.R
  
  data(Duncan, package="car")
  
  mod.1 <- lm(prestige ~ type + poly(income, degree=2, raw=TRUE), data=Duncan)
  

  # (2) focal: polynomial, constant: factor
  
  print(plot(Effect(c("income"), mod.1, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("income"), mod.1, partial.residual=TRUE)$fit, 
                        Effect(c("income"), mod.1, xlevels=100)$fit)))
    stop("failed test 2 (2)")
  
  
  # (3) focal: factor*polynomial, constant: polynomial
  
  mod.2 <- lm(prestige ~ type*poly(income, degree=2, raw=TRUE) + 
                poly(education, degree=2, raw=TRUE), data=Duncan)

  print(plot(Effect(c("type", "income"), mod.2, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("type", "income"), mod.2, partial.residual=TRUE)$fit, 
                        Effect(c("type", "income"), mod.2, xlevels=100)$fit)))
    stop("failed test 2 (3)")
  
  
  # (4) focal: polynomial, constant: factor*polynomial
  
  print(plot(Effect(c("education"), mod.2, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("education"), mod.2, partial.residual=TRUE)$fit, 
                        Effect(c("education"), mod.2, xlevels=100)$fit)))
    stop("failed test 2 (4)")  
  
  
  # (6) focal: factor*polynomial, constant: polynomial, using predict() & orthog. polys.
  
  mod.3 <- lm(prestige ~ type*poly(income, degree=2) + poly(education, degree=2), data=Duncan)
  
  print(plot(Effect(c("type", "income"), mod.3, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("type", "income"), mod.3, partial.residual=TRUE)$fit, 
                        Effect(c("type", "income"), mod.3, xlevels=100)$fit)))
    stop("failed test 2 (6)")  
  
  # (7) focal: factor, constant: poly*poly
  
  mod.4 <- lm(prestige ~ type + poly(income, 2)*poly(education, 2), data=Duncan)
  
  print(plot(Effect(c("income", "education"), mod.4, partial.residuals=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("income", "education"), mod.4, partial.residuals=TRUE)$fit, 
                        Effect(c("income", "education"), mod.4, 
                               xlevels=list(income=100, 
                                            education=quantile(Duncan$education, probs=seq(0.2, 0.8, by=0.2))))$fit)))
    stop("failed test 2 (7)")  
  
  # (9) focal: covariate, constant: 2 factors and 1 covariate, 3-way interaction
  
  data(Mroz, package="car")
  mod.6 <- lm(lwg ~ inc + age*hc*wc, data=Mroz)

  print(plot(Effect(c("inc"), mod.6, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("inc"), mod.6, partial.residual=TRUE)$fit, 
                        Effect(c("inc"), mod.6, xlevels=100)$fit)))
    stop("failed test 2 (9-1)")  
  
  print(plot(Effect(c("age", "hc", "wc"), mod.6, partial.residual=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("age", "hc", "wc"), mod.6, partial.residual=TRUE)$fit, 
                        Effect(c("age", "hc", "wc"), mod.6, xlevels=100)$fit)))
    stop("failed test 2 (9-2)")  
  
  # additional tests of partial residuals
  
  mod.7 <- lm(prestige ~ income*type + education, data=Prestige)
  print(plot(Effect(c("income", "type"), mod.7, partial.residuals=TRUE), show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("income", "type"), mod.7, partial.residuals=TRUE)$fit, 
                        Effect(c("income", "type"), mod.7, xlevels=100)$fit)))
    stop("failed test 2 (additional-1)")  
  
  Mroz2 <- Mroz
  Mroz2$hc <- as.numeric(Mroz$hc) - 1
  Mroz2$wc <- as.numeric(Mroz$wc) - 1
  mod.8 <- lm(lwg ~ inc*age*k5 + hc*wc, data=Mroz2)
  print(plot(Effect(c("inc", "age", "k5"), mod.8, partial.residuals=TRUE, xlevels=list(k5=0:1)), 
             show.fitted=TRUE))
  
  if (!isTRUE(all.equal(Effect(c("inc", "age", "k5"), mod.8, partial.residuals=TRUE, xlevels=list(k5=0:1))$fit, 
                        Effect(c("inc", "age", "k5"), 
                               mod.8, partial.residuals=TRUE, 
                               xlevels=list(k5=0:1, inc=100, age=quantile(Mroz2$age, seq(.2, .8, by=.2))))$fit)))
    stop("failed test 2 (additional-2)")  
  
  print(plot(Effect(c("hc", "wc"), mod.8, partial.residuals=TRUE, xlevels=list(hc=0:1, wc=0:1)), 
       show.fitted=TRUE, smooth.residuals=FALSE, residuals.pch="."))
}


