
 library(rcdd)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)

 out <- scdd(qux)
 print(out)

 out <- scdd(out$output)
 print(out)

 out <- scdd(out$output, roworder = "mincutoff")
 print(out)

 ##########################

 quux <- addHeq(1:d, (d + 1) / 2, qux)
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

 ##########################

 quux <- makeH(a2 = rep(1, d), b2 = 1)
 print(quux)

 out <- scdd(quux)
 print(out)

 ##########################

 quux <- qux[- 1, ]
 print(quux)

 out <- scdd(quux)
 print(out)

 ##########################

 quux[ , 2] <- quux[ , 2] + 1
 print(quux)

 print(quux)
 out <- scdd(quux)
 print(out)

