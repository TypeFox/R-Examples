require(lmerTest)

carrots$Income <- as.factor(carrots$Income)
carrots$Age <- as.factor(carrots$Age)

m.carrots <- lmer(Preference ~ Age*Homesize*Income
                  +(1+sens2|Consumer), data=carrots)

an.carrots <- anova(m.carrots) ## OK with SAS

TOL <- 1e-4 # for the check

#with 4 decimals should agree with SAS output
#numbers vefore decimals should agree with SAS output
stopifnot(
  all.equal(an.carrots[,"Pr(>F)"], c(0.3541, 0.0838, 0.2879, 0.3263, 0.6173, 0.2652, NA), tol = TOL), 
  all.equal(round(an.carrots$DenDF), c(77, 77, 77, 77, 77, 77, 0))
  , TRUE)


tools::assertError(step(m.carrots)) ## error because of the convergence issues

load(system.file("testdata","bread.RData", package="lmerTest"))
bread$Children <- as.factor(bread$Children)
bread$Exposure <- as.factor(bread$Exposure)
bread$type <- as.factor(bread$type)
bread$extra <- as.factor(bread$extra)

m.bread <- lmer(Liking ~ ref + type * extra + (1 | Exposure) + (1 | Children), data=bread)
an.bread <- anova(m.bread) ## OK with SAS

stopifnot(
  all.equal(an.bread[,"Pr(>F)"], c(NA, 1e-7, 0.5864, 0.2783), tol = TOL), 
  all.equal(round(an.bread$DenDF), c(0, 1798, 1796, 1796))
  , TRUE)

s.bread <- step(m.bread) ## returns lsmeans and difflsmeans, in SAS -non-est, why?!! compare with lsmeans/ doby packages


# results for lsmeans for type 1-2 should agree with SAS
stopifnot(all.equal(s.bread$anova.table[, "F.value"],
                    c(1.1762, 0.2983, NA, 25.1939), tolerance=TOL, check.attributes=FALSE),  
          all.equal(s.bread$diffs.lsmeans.table["type 1 - 2",], 
                    c(0.3599, 0.0717, 1799.7, 5.02, 0.2193, 0.5006, 0),      
                    tolerance=TOL, check.attributes=FALSE),                     
                    TRUE)

