# Source file that executes the tests that have to run only once
require(testthat)
require(randomizeR)
# path <- "U:/IDeAl/IDeAlGrit/randomizeR/randomizeR/tests/testthat/OnlyOnce/"
# path <- "D:\\David_lokal\\Github\\Rpackage\\randomizeR\\randomizeR\\tests\\testthat\\OnlyOnce\\"
# path <- "C:/Users/duschner/package/randomizeR/randomizeR/tests/testthat/"
setwd(path)

liste <- Sys.glob("./OnlyOnce/*.R")
for(i in liste) {
  test_file(path = i)
}
