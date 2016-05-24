###
### $Id: midwt.R 25 2014-06-20 21:02:46Z plroebuck $
###

options(warn=1)
library(rwt)


##-----------------------------------------------------------------------------
test.midwt <- function(input, expected) {
   ret.mdwt <- rwt::mdwt(input$signal, input$filter, input$nlevels)
   y <- ret.mdwt$y
   L <- ret.mdwt$L
   result <- rwt::midwt(y, input$filter, L)
   identical(all.equal(result,
                       expected,
                       tolerance=0.000001),
             TRUE)
}


sig <- rwt::makesig(SIGNAL.LIN.CHIRP, 8)
h <- rwt::daubcqf(4, PHASE.MINIMUM)
midwt.expected <- list(x = matrix(data = c( 0.04906767,
                                            0.1950903,
                                            0.4275551,
                                            0.7071068,
                                            0.941544,
                                            0.9807853,
                                            0.671559,
                                            0.0000),
                                 nrow = 1),
                      L = 1)

test.midwt(list(signal = sig$x, filter = h$h.0, nlevels = 1), midwt.expected)

