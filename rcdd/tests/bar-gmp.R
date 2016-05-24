
 library(rcdd)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
 qux <- d2q(qux)
 print(qux)

 out <- scdd(qux)
 print(out)

 out <- scdd(out$output)
 print(out)

