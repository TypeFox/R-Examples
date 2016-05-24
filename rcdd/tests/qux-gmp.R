
 library(rcdd)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), "1")
 print(qux)

 out <- scdd(qux)
 print(out)

 ##########################

 quux <- addHeq(1:d, z2q(d + 1, 2), qux)
 print(quux)

 out <- scdd(quux)
 print(out)

 ##########################

 quux <- addHin(c(1, -1, 0, 0), 0, qux)
 print(quux)

 out <- scdd(quux)
 print(out)

 ##########################

 quux <- qux[- 4, ]
 print(quux)

 out <- scdd(quux)
 print(out)

