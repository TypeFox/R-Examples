require(lmerTest)

testType1 <- TRUE
##
## compare output with SAS for rank deficient models
##
load(system.file("testdata","cltlike.RData", package="lmerTest"))

if(testType1){

fm2 <- lmer(liking ~ pig.type + group + landro.male + lskatol.male +
              sens1 + sens2 + gender + agegroup +
              sens.ska + sens.and +
              sens.ska:lskatol + sens.and:landro +
              (1|pig) + (1|consumer), data=cltlike)

## for type3 hypothesis matrix
anType3 <- anova(fm2)

TOL <- 1e-3 # for the check

#with 3 decimals should agree with SAS output
#numbers vefore decimals should agree with SAS output
stopifnot(
  all.equal(anType3[,"Pr(>F)"], c(NA, 0.28047, NA, NA, 0.14458, 0.06193, 0.31344,
                                  0.07603, 0.14608, 0.07304, 0.90075, 0.04722), tol = TOL), 
  all.equal(round(anType3$DenDF), c(0, 55, 0, 0, 61, 66, 170, 169, 477, 192, 709, 692))
  , TRUE)

## for type 1 hypothesis matrix
anType1 <- anova(fm2, method.grad="Richardson", type=1)

#with 3 decimals should agree with SAS output
#numbers vefore decimals should agree with SAS output
stopifnot(
  all.equal(anType1[,"Pr(>F)"], c(0.012235, 0.286167, 0.767076, 0.005902, 0.400181,
                                  0.05243, 0.10677, 0.106436, 0.034507, 0.14896, 
                                  0.842839, 0.047221), tol = TOL), 
  all.equal(round(anType1$DenDF), c(62, 53, 54, 58, 57, 66, 169, 169, 172, 168, 710, 692))
  , TRUE)


##
## compare output with SAS for rank deficient models
##
}