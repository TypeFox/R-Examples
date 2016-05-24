
 library(rcdd)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)
 print(qux)
 out <- scdd(qux)
 print(out)

 qux1 <- addHeq(1:d, 2.2, qux)
 print(qux1)
 out1 <- scdd(qux1)
 print(out1)

 qux1q <- qux1
 qux1q[6, 2] <- "11/5"
 print(qux1q)
 out1q <- scdd(qux1q)
 print(out1q)

 out1$output
 q2d(out1q$output)
 qmq(out1$output, out1q$output)

