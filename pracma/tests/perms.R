##
##  p e r m s . R  Test Suite
##

perms <- pracma::perms
randperm <- pracma::randperm

identical(perms(2), matrix(2, 1, 1))
identical(perms(c(1, 2)), matrix(c(2, 1, 1, 2), 2, 2))
identical(nrow(perms(1:6)), 720L)

all(c(1,2,3) %in% randperm(c(1,2,3)))
