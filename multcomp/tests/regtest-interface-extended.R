### test the extended interpreter for the left hand side of linear hypotheses
### Features: 
### - fully recursive expression parser built upon a small code fragment copied over from base::codetools
### - the parser stops if any of the following conditions is not met:
###   - any variable  must be addressed only once
###   - all operators and functions must finally evaluate to a a real valued literal
###   - function parameters must not denote an effect name
###   - effects can not be multiplied or divided by another effect
###   - additive or subtractive terms involving an effect and a numeric 
###     constants must not be specified 
###   - coefficients associated with named effects must not evaluate to zero 
###
### Examples:
###   x1 + x1             == 0  -> not accepted
###   x1 + x2 -1          == 0  -> not accepted
###   x1 * x2             == 0  -> not accepted
###   x1 / x2             == 0  -> not accepted
###   f(x1)               == 0  -> not accepted if x1 denotes an effect
###   2*3                 == 6  -> not accepted because no effect was named
###   x1 + x2*0           == 0  -> not accepted because this is likely an oversight
###   x1 + 3*(4-5+1)*x2   == 0  -> not accepted because this is likely an oversight
###   x1*3/0              == 0  -> not accepted because coefficient would become infinite
###   log(-1)*x1          == 0  -> not accepted, because the result is not finite
###   x1 + x2 +0          == 0  -> accepted because adding zero does not make a difference
###   sin(pi/2) * x1      == 0  -> accepted if 'pi' is not an effect
###   sin(Pi/2) * x1      == 0  -> accepted if 'Pi' is not an effect. However, if the environment does not define Pi the evaluation may still fail.


tmp <- multcomp:::chrlinfct2matrix( c( l01 = " x1 - x2 = 2"
                                     , l02 = " x2 + 3 * x3 = 1"
                                     , l03 = " (x1 - x2) - (x3 - x4) =  0"
                                     , l04 = "+(x1 - x2)*-2 - (1/3+2)*( +x3 - 2*x4 ) = -1" 
                                     , l05 = "-(x1 - x2)*-2 - (1/3+2)*( -x3 - 2*x4 ) = -2" 
                                     , l06 = "-(x1 - x2)*-2  - (1/3+2)*( -x3 - 2*x4 )*7/-10 = -3" 
                                     , l07 = "-1*(x1:x2 - x1:x2:x3) - x3 = -4"
                                     , l08 = "-(x1:x2 - x1:x2:x3) - x3 = -4"
                                     , l09 = "-(x1:x2 - 3*x1:x2:x3)*-2 - x3 -5/3*-x4= -5"
                                     , l10 = "--cos(pi/2)*x1 - 10*(log(10^-3)+1)*-x2 -10^-3*x3 + -exp(-2)*x4= -6"
                                     , l11 = " x1 + x2 + 0 = -7"
                                     ),  c('x1','x2','x3','x4','x1:x2','x1:x2:x3') )

stopifnot(max(abs( dK <- tmp$K - 
                         rbind( c(           1,                 -1,                0,                 0,       0,     0 )
                              , c(           0,                  1,                3,                 0,       0,     0 )
                              , c(           1,                 -1,               -1,                 1,       0,     0 )
                              , c(          -2,                  2,         -(1/3+2),         2*(1/3+2),       0,     0 )
                              , c(           2,                 -2,          (1/3+2),         2*(1/3+2),       0,     0 )
                              , c(           2,                 -2,    (1/3+2)*-7/10,   2*(1/3+2)*-7/10,       0,     0 )
                              , c(           0,                  0,               -1,                 0,      -1,     1 )
                              , c(           0,                  0,               -1,                 0,      -1,     1 )
                              , c(           0,                  0,               -1,           -5/3*-1,       2,    -6 )
                              , c( --cos(pi/2),  10*(log(10^-3)+1),           -10^-3,           -exp(-2),      0,     0 )
                              , c(           1,                  1,                0,                 0,       0,     0 )
                              ))) < sqrt(.Machine$double.eps))

stopifnot(max(abs( tmp$m - 
                   c(  2
                    ,  1
                    ,  0
                    , -1
                    , -2
                    , -3
                    , -4
                    , -4
                    , -5
                    , -6
                    , -7
                    ))) < sqrt(.Machine$double.eps))

expectFail <- function(testname, x) {
 if ( class(x) != 'try-error' ) {
      stop(testname, ' unexpectedly succeeded. Result is: ', paste(x, collapse = ', '),'\n')
 }
 message(testname, ' expectedly failed. Message is: ', attr(x,'condition')$message, '\n')
}

expectFail('test 01',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 - x1  = 0"), c('x1','x2')), silent=T))

expectFail('test 02',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 - X2  = 0"), c('x1','x2')), silent=T))

expectFail('test 03',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 - x2 -1 = 0"), c('x1','x2')), silent=T))

expectFail('test 04',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 * x2  = 0"), c('x1','x2')), silent=T))

expectFail('test 05',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 / x2  = 0"), c('x1','x2')), silent=T))

expectFail('test 06',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 - exp(x2)  = 0"), c('x1','x2')), silent=T))

expectFail('test 07',  try( multcomp:::chrlinfct2matrix(c(l1 = "sin(Pi)*x1   = 0"), c('x1','x2')), silent=T))

expectFail('test 08',  try( multcomp:::chrlinfct2matrix(c(l1 = "3*4 = 0"), c('x1','x2')), silent=T))

expectFail('test 09',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1 + 3*(4-5+1)*x2 = 0"), c('x1','x2')), silent=T))

expectFail('test 10',  try( multcomp:::chrlinfct2matrix(c(l1 = "x1*3/0 = 0"), c('x1','x2')), silent=T))

expectFail('test 11',  try( multcomp:::chrlinfct2matrix(c(l1 = "log(-1)*x1 = 0"), c('x1','x2')), silent=T))
