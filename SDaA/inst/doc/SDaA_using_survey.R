### R code from vignette source 'SDaA_using_survey.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: config
###################################################
	options(width=65)


###################################################
### code chunk number 2: loadPkg
###################################################
  library(SDaA)
  library(survey)


###################################################
### code chunk number 3: SDaA_using_survey.Rnw:45-49
###################################################
  ### Example 3.2, p. 63 
  agsrsDesign <- svydesign(ids=~1, weights = ~1, data = agsrs)
  svyratio(numerator = ~acres92, denominator = ~acres87, 
      design = agsrsDesign) # proportion B hat


###################################################
### code chunk number 4: acreagetwoyears
###################################################
  ### part of Example 3.2, p. 64
  plot(I(acres92/10^6) ~ I(acres87/10^6), 
      xlab = "Millions of Acres Devoted to Farms (1987)", 
      ylab = "Millions of Acres Devoted to Farms (1992)", data = agsrs)
  abline(lm(I(acres92/10^6) ~ 0 + I(acres87/10^6), # through the origin 
          data = agsrs), col = "red", lwd = 2)


###################################################
### code chunk number 5: seedlingData
###################################################
  ### Example 3.5, p. 72, table 3.3
  seedlings <- data.frame(tree = 1:10, 
      x = c(1, 0, 8, 2, 76, 60, 25, 2, 1, 31),
      y = c(0, 0, 1, 2, 10, 15, 3, 2, 1, 27))
  names(seedlings) <- c("tree", "x", "y")


###################################################
### code chunk number 6: seedlingsplot
###################################################
  plot(y ~ x, data = seedlings, xlab = "Seedlings Alive (March 1992)",
      ylab = "Seedlings That Survived (February 1994)")
  # abline(lm(y ~ 0 + x, data = seedlings), lwd = 2, col = "red")
  # TODO: add proper abline


###################################################
### code chunk number 7: SDaA_using_survey.Rnw:93-99
###################################################
  ### Example 3.6, p. 75
  pf <- data.frame(photo = c(10, 12, 7, 13, 13, 6, 17, 
          16, 15, 10, 14, 12, 10, 5,
          12, 10, 10, 9, 6, 11, 7, 9, 11, 10, 10),
      field = c(15, 14, 9, 14, 8, 5, 18, 15, 13, 15, 11, 15, 12,
          8, 13, 9, 11, 12, 9, 12, 13, 11, 10, 9, 8))


###################################################
### code chunk number 8: modelreg
###################################################
	### Example 3.9, p. 83
  recacr87 <- agsrs$acres87
  recacr87[recacr87 > 0] <- 1/recacr87[recacr87 > 0] # cf. p. 450
  model1 <- lm(acres92 ~ 0 + acres87, weights = recacr87, data = agsrs)
  summary(model1)


###################################################
### code chunk number 9: modelregplot
###################################################
  ### Figure 3.6, p. 85
  wtresid <- resid(model1) / sqrt(agsrs$acres87) 
  plot(wtresid ~ I(agsrs$acres87/10^6), 
      xlab = "Millions of Acres Devoted to Farms (1987)",
      ylab = "Weighted Residuals")


###################################################
### code chunk number 10: acresboxplot
###################################################
	boxplot(acres92/10^6 ~ region, xlab = "Region", 
      ylab = "Millions of Acres", data = agstrat)


###################################################
### code chunk number 11: gpaex
###################################################
  ### Example 5.2, p. 137 middle
  GPA <- cbind(expand.grid(1:4, 1:5), 
      gpa = c(3.08, 2.60, 3.44, 3.04, 2.36, 3.04, 3.28, 2.68, 2.00, 2.56, 
          2.52, 1.88, 3.00, 2.88, 3.44, 3.64, 2.68, 1.92, 3.28, 3.20))
  names(GPA)[1:2] <- c("person_num", "cluster")
  GPA$pwt <- 100/5
  
  clusterDesign <- svydesign(ids = ~ cluster, weights = ~ pwt, data = GPA)
  
  svytotal(~ gpa, design = clusterDesign)
  #      total     SE 
  # gpa 1130.4 67.167
  
  # Stata results: 1130.4   67.16666 ---> corresponds perfectly


###################################################
### code chunk number 12: SDaA_using_survey.Rnw:172-175
###################################################
  ### Figure 5.3
	plot(volume ~ clutch, xlim = c(0,200), pch=19, data = coots,
      xlab = "Clutch Number", ylab = "Egg Volume")


###################################################
### code chunk number 13: SDaA_using_survey.Rnw:178-181
###################################################
  ### Figure 5.3
  plot(volume ~ clutch, xlim = c(0,200), pch=19, data = coots,
      xlab = "Clutch Number", ylab = "Egg Volume")


###################################################
### code chunk number 14: readStatepop
###################################################
  data(statepop)
  statepop$psi <- statepop$popn / 255077536


###################################################
### code chunk number 15: SDaA_using_survey.Rnw:197-201
###################################################
	### page 191, figure 6.1
  plot(phys ~ psi, data = statepop, 
       xlab = expression(paste(Psi[i], " for County")),
       ylab = "Physicians in County (in thousands)")


###################################################
### code chunk number 16: SDaA_using_survey.Rnw:210-215
###################################################
	### Figure 7.1
  data(htpop)
  popecdf <- ecdf(htpop$height)
  plot(popecdf, do.points = FALSE, ylab = "F(y)", 
       xlab = "Height Value, y")


###################################################
### code chunk number 17: SDaA_using_survey.Rnw:218-223
###################################################
  ### Figure 7.2
  minht <- min(htpop$height)
  breaks <- c(minht-1, seq(from = minht, to = max(htpop$height), by = 1))
  hist(htpop$height, ylab = "f(y)", breaks = breaks, 
       xlab = "Height Value, y", freq = FALSE)


###################################################
### code chunk number 18: SDaA_using_survey.Rnw:226-230
###################################################
  ### Figure 7.3
  data(htsrs)
  hist(htsrs$height, ylab = "Relative Frequency", 
       xlab = "Height (cm)", freq = FALSE)


###################################################
### code chunk number 19: SDaA_using_survey.Rnw:235-239
###################################################
  ### Figure 7.4
  data(htstrat)
  hist(htstrat$height, ylab = "Relative Frequency", 
      xlab = "Height (cm)", freq = FALSE)


###################################################
### code chunk number 20: SDaA_using_survey.Rnw:242-247
###################################################
	### Figure 7.5 (a)
  minht <- min(htstrat$height)
  breaks <- c(minht-1, seq(from = minht, to = max(htstrat$height), by = 1))
  hist(htstrat$height, ylab = expression(hat(f)(y)), breaks = breaks, 
       xlab = "Height Value, y", freq = FALSE)


###################################################
### code chunk number 21: SDaA_using_survey.Rnw:250-254
###################################################
  ### Figure 7.5 (b)
  stratecdf <- ecdf(htstrat$height)
  plot(stratecdf, do.points = FALSE, ylab = expression(hat(F)(y)), 
       xlab = "Height Value, y")


###################################################
### code chunk number 22: SDaA_using_survey.Rnw:259-262
###################################################
	### Figure 7.6
  data(syc)
  hist(syc$age, freq = FALSE, xlab = "Age")


###################################################
### code chunk number 23: SDaA_using_survey.Rnw:273-286
###################################################
  ### Figure 7.8
  sycdesign <- svydesign(ids= ~ psu, strata = ~ stratum,
     data = syc, weights=~finalwt)
  # p. 235: "Each of the 11 facilities with 360 or more youth
  # formed its own stratum (strata 6-16)", so in order
  # to avoid a lonely.psu error message
  #  Error in switch(lonely.psu, certainty = scale * crossprod(x), remove = scale *  : 
  #          Stratum (6) has only one PSU at stage 1
  # we set the option to "certainty" for this example
  # to see the problem, use: by(syc$psu, syc$stratum, unique) 
  oo <- options(survey.lonely.psu = "certainty")
  svyboxplot(age ~ factor(stratum), design = sycdesign) # mind the factor
  options(oo)


###################################################
### code chunk number 24: SDaA_using_survey.Rnw:293-298
###################################################
  ### Figure 7.9
  library(ggplot2)
  p <- ggplot(syc, aes(x = factor(stratum), y = factor(age)))
  g <- p + stat_sum(aes(group=1, weight = finalwt, size = ..n..)) 
  print(g)


###################################################
### code chunk number 25: SDaA_using_survey.Rnw:310-315
###################################################
	### Figure 7.10
  oo <- options(survey.lonely.psu = "certainty")
  sycstrat5 <- subset(sycdesign, stratum == 5)
  svyboxplot(age ~ factor(psu), design = sycstrat5)
  options(oo)


###################################################
### code chunk number 26: SDaA_using_survey.Rnw:318-323
###################################################
  ### Figure 7.11
  sycstrat5df <- subset(syc, stratum == 5)
  p <- ggplot(sycstrat5df, aes(x = factor(psu), y = factor(age)))
  g <- p + stat_sum(aes(group=1, weight = finalwt, size = ..n..)) 
  print(g)


###################################################
### code chunk number 27: SDaA_using_survey.Rnw:341-348
###################################################
  ### Example 10.1
  hh <- rbind(c(119, 188),
              c(88, 105))
  rownames(hh) <- c("cableYes", "cableNo")
  colnames(hh) <- c("computerYes", "computerNo")
  addmargins(hh)
  chisq.test(hh, correct = FALSE)  # OK      


###################################################
### code chunk number 28: SDaA_using_survey.Rnw:354-364
###################################################
	### Example 10.2 (nursing students and tutors)
  nst <- rbind(c(46, 222),
               c(41, 109),
               c(17, 40),
               c(8, 26))
  colnames(nst) <- c("NR", "R")
  rownames(nst) <- c("generalStudent", "generalTutor", "psychiatricStudent", 
      "psychiatricTutor")
  addmargins(nst)
  chisq.test(nst, correct = FALSE) # OK        


###################################################
### code chunk number 29: SDaA_using_survey.Rnw:367-376
###################################################
	### Example 10.3 (Air Force Pilots)
  afp <- data.frame(nAccidents = 0:7, 
                    nPilots = c(12475, 4117, 1016, 269, 53, 14, 6, 2))
  # estimate lambda
  lambdahat <- sum(afp$nAccidents * afp$nPilots  / sum(afp$nPilots))
  # expected counts
  observed <- afp$nPilots
  expected <- dpois(0:7, lambda = lambdahat) * sum(afp$nPilots)
  sum((observed - expected)^2 / expected) # NOT OK


###################################################
### code chunk number 30: SDaA_using_survey.Rnw:383-390
###################################################
	### Example 10.4
  hh2 <- rbind(c(238, 376),
               c(176, 210))
  rownames(hh2) <- c("cableYes", "cableNo")
  colnames(hh2) <- c("computerYes", "computerNo")
  addmargins(hh2)
  chisq.test(hh2, correct = FALSE)  # OK


###################################################
### code chunk number 31: SDaA_using_survey.Rnw:399-400
###################################################
	### example 10.5


