## follow http://adv-r.had.co.nz/Testing.html
print(sprintf("Test seqminer [version %s]", packageVersion("seqminer")))

print("Platform")
print(str(.Platform))

print("Sys.info()")
print(Sys.info())

print(citation("seqminer"))

library(testthat)
library(seqminer)
## test code are under inst/tests
## test_package("seqminer", reporter="tap")
test_package("seqminer")
