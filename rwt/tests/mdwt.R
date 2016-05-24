###
### $Id: mdwt.R 25 2014-06-20 21:02:46Z plroebuck $
###

options(warn=1)
library(rwt)


##-----------------------------------------------------------------------------
test.mdwt <- function(input, expected) {
   result <- rwt::mdwt(input$signal, input$filter, input$nlevels)
   identical(all.equal(result,
                       expected,
                       tolerance=0.000001),
             TRUE)
}


sig <- rwt::makesig(SIGNAL.LIN.CHIRP, 8)
h <- rwt::daubcqf(4, PHASE.MINIMUM)
mdwt.expected <- list(y = matrix(data = c( 1.109692,
                                           0.8766618,
                                           0.8203919,
                                          -0.5200741,
                                          -0.03392767,
                                           0.1001107,
                                           0.2200882,
                                          -0.1400816),
                                 nrow = 1),
                      L = 2)

test.mdwt(list(signal = sig$x, filter = h$h.0, nlevels = 2), mdwt.expected)

