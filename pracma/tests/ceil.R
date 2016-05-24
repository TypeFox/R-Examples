##
##  c e i l . R  tests
##

ceil <- pracma::ceil
Fix <- pracma::Fix

identical(ceil(0), 0)
identical(ceil(-1), -1)
identical(ceil(-1.5), -1)
identical(ceil(1), 1)
identical(ceil(1.5), 2)

identical(Fix(0), 0)
identical(Fix(-1), -1)
identical(Fix(-1.5), -1)
identical(Fix(1), 1)
identical(Fix(1.5), 1)
