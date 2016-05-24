require(testthat)
require(attribrisk)

data(stroke.dat)
smoke.baseline <- stroke.dat$smoke
smoke.baseline[smoke.baseline != "Never"] <- "Never"

DBP.baseline <- .90 * stroke.dat$diastolic
stroke.target1 <- data.frame(diastolic = DBP.baseline, smoke = smoke.baseline)

set.seed(1013)
x <- attribrisk(cases ~ age + expos(smoke) + expos(diastolic), 
                data = stroke.dat, varmethod='jackknife',
                baseline=stroke.target1)

expect_equal(x$attribrisk,0.4282319, 
             label=paste("Attribrisk estimate not within tolerance for following call: ", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


expect_equal(object=sqrt(x$var), 0.040287944, 
             label=paste("Attribrisk estimate for sterr not in tolerance for following call:", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


