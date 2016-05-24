
 library(rcdd)

 d <- 4
 qux <- makeH(- diag(d), rep(0, d), rep(1, d), 1)

 out <- scdd(qux, adjacency = TRUE)
 print(out)

 ##########################

 quux <- addHeq(1:d, (d + 1) / 2, qux)
 print(quux)

 out <- scdd(quux, adjacency = TRUE)
 print(out)

 out <- scdd(quux, incidence = TRUE)
 print(out)

 out <- scdd(quux, inputadjacency = TRUE)
 print(out)

 out <- scdd(quux, inputincidence = TRUE)
 print(out)

