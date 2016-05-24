require(testthat)
require(attribrisk)

data(chapter.dat)


set.seed(1010)
x <- attribrisk(cases ~ expos(hbp), data = chapter.dat, varmethod='jackknife')

expect_equal(x$attribrisk, 0.3130583, 
             label=paste("Attribrisk estimate not within tolerance for following call: ", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


expect_equal(object=sqrt(x$var),  0.03834875, 
             label=paste("Attribrisk estimate for sterr not in tolerance for following call:", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)

