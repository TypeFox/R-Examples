##
##  m o d . R  tests
##

mod <- pracma::mod
rem <- pracma::rem

identical(mod(0, 0), 0)
identical(mod(1, 0), 1)
identical(mod(0, 2), 0)
identical(mod(5, 3), 2)
identical(mod(5, -3), -1)
identical(mod(-5, 3), 1)
identical(mod(-5, -3), -2)

identical(rem(0, 0), NaN)
identical(rem(1, 0), NaN)
identical(rem(0, 2), 0)
identical(rem(5, 3), 2)
identical(rem(5, -3), 2)
identical(rem(-5, 3), -2)
identical(rem(-5, -3), -2)
