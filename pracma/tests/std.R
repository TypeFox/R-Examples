##
##  s t d . r  tests
##

std <- pracma::std
std_err <- pracma::std_err

all.equal(std(1:10), 3.0277, tolerance=0.0001)
all.equal(std(1:10, flag=0), 3.0277, tolerance=0.0001)
all.equal(std(1:10, flag=1), 2.8723, tolerance=0.0001)

all.equal(std_err(1:10), 0.9574271, tolerance=0.0001)
