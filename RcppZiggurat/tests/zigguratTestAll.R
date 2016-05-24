
library(RcppZiggurat)

cat("Default\n")
zsetseed(123456890)
print(zrnorm(10), digits=12)

cat("MT\n")
zsetseedMT(123456890)
print(zrnormMT(10), digits=12)

cat("LZLLV\n")
zsetseedLZLLV(123456890)
print(zrnormLZLLV(10), digits=12)

cat("GSL\n")
zsetseedGSL(123456890)
print(zrnormGSL(10), digits=12)

cat("QL\n")
zsetseedQL(123456890)
print(zrnormQL(10), digits=12)

cat("Gretl\n")
zsetseedGl(123456890)
print(zrnormGl(10), digits=12)

cat("R\n")
set.seed(1234567890)
print(zrnormR(10), digits=12)
