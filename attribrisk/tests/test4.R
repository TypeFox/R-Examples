require(testthat)
require(attribrisk)

data(stroke.dat)


smoke.baseline <- stroke.dat$smoke
smoke.baseline[smoke.baseline == "Current"] <- "Former"
DBP.baseline <- ifelse(stroke.dat$diastolic > 150,
                       .75 * stroke.dat$diastolic, stroke.dat$diastolic)
DBP.baseline <- ifelse(stroke.dat$diastolic > 120 & stroke.dat$diastolic <= 150,
                       .85 * stroke.dat$diastolic, DBP.baseline)
DBP.baseline <- ifelse(stroke.dat$diastolic > 100 & stroke.dat$diastolic <= 120,
                       .90 * stroke.dat$diastolic, DBP.baseline)
DBP.baseline <- ifelse(stroke.dat$diastolic > 85 & stroke.dat$diastolic <= 100,
                       .95 * stroke.dat$diastolic, DBP.baseline)
stroke.target2 <- data.frame(diastolic = DBP.baseline, smoke = smoke.baseline)

set.seed(1014)
x <- attribrisk(cases ~ age + expos(smoke) + expos(diastolic), 
                data = stroke.dat, baseline = stroke.target2, 
                varmethod = 'jackknife')

expect_equal(x$attribrisk,0.1170187, 
             label=paste("Attribrisk estimate not within tolerance for following call: ", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


expect_equal(object=sqrt(x$var), 0.01829625, 
             label=paste("Attribrisk estimate for sterr not in tolerance for following call:", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)

