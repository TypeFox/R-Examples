
 library(rcdd)

 # needed because now uses R RNG for matrix row permutation
 set.seed(42)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
 print(qux)

 out <- scdd(qux)
 print(out)

 out <- scdd(out$output)
 print(out)

 out <- scdd(qux, roworder = "randomrow")
 print(out)

 out <- scdd(qux, roworder = "maxcutoff")
 print(out)

