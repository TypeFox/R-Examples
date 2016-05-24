###
### $Id: mirdwt.R 25 2014-06-20 21:02:46Z plroebuck $
###

options(warn=1)
library(rwt)


##-----------------------------------------------------------------------------
test.mirdwt <- function(input, expected) {
   ret.mrdwt <- rwt::mrdwt(input$signal, input$filter, input$nlevels)
   yl <- ret.mrdwt$yl
   yh <- ret.mrdwt$yh
   L  <- ret.mrdwt$L
   result <- rwt::mirdwt(yl, yh, input$filter, L)
   identical(all.equal(result,
                       expected,
                       tolerance=1e-17),
             TRUE)
}


sig <- rwt::makesig(SIGNAL.LEOPOLD, 8)
h <- rwt::daubcqf(4, PHASE.MINIMUM)
mirdwt.expected <- list(x = matrix(data = c( 8.6736e-018,
                                             1.0000,
                                             8.6736e-018,
                                            -1.3878e-017,
                                             0,
                                             0,
                                             0,
                                            -1.3878e-017),
                                 nrow = 1),
                      L = 1)

test.mirdwt(list(signal = sig$x, filter = h$h.0, nlevels = 1), mirdwt.expected)

