##
##  f i n d i n t e r v a l s . R  Test suite
##

findintervals <- pracma::findintervals

identical(findintervals(0, zapsmall(sin(seq(0, 10*pi, len=100)))),
          as.integer(c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)))
