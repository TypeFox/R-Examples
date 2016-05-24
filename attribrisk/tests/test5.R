require(testthat)
require(attribrisk)

data(benichou)

set.seed(1015)
x <- attribrisk(cases ~ expos(alcohol80), data = benichou, 
                 varmethod='jackknife')

expect_equal(x$attribrisk,0.3948949, 
             label=paste("Attribrisk estimate not within tolerance for following call: ", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)


expect_equal(object=sqrt(x$var),  0.05028568, 
             label=paste("Attribrisk estimate for sterr not in tolerance for following call:", paste(deparse(x$call), collapse="")), tolerance=1.0e-6, scale=1)
