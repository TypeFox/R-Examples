###
### $Id: mrdwt.R 25 2014-06-20 21:02:46Z plroebuck $
###

options(warn=1)
library(rwt)


##-----------------------------------------------------------------------------
test.mrdwt <- function(input, expected) {
   result <- rwt::mrdwt(input$signal, input$filter, input$nlevels)
   identical(all.equal(result,
                       expected,
                       tolerance=0.0000001),
             TRUE)
}


sig <- rwt::makesig(SIGNAL.LEOPOLD, 8)
h <- rwt::daubcqf(4, PHASE.MINIMUM)
mrdwt.expected <- list(yl = matrix(data = c( 0.8365163,
                                             0.4829629,
                                             0,
                                             0,
                                             0,
                                             0,
                                            -0.1294095,
                                             0.2241439)),
                       yh = matrix(data = c(-0.2241439,
                                            -0.1294095,
                                             0,
                                             0,
                                             0,
                                             0,
                                            -0.4829629,
                                             0.8365163)),
                       L = 1)

test.mrdwt(list(signal = sig$x, filter = h$h.0, nlevels = 1), mrdwt.expected)

