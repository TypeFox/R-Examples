require(testthat)
require(attribrisk)

data(benichou)

set.seed(1016)
x <- attribrisk(cases ~ expos(alcohol), data = benichou, 
                 varmethod='jackknife')

expect_equal(x$attribrisk,0.7088729, 
             label=paste("Attribrisk estimate not within tolerance for following call: ", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


expect_equal(object=sqrt(x$var), 0.05454432, 
             label=paste("Attribrisk estimate for sterr not in tolerance for following call:", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)
