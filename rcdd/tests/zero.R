
 library(rcdd)

 qux <- matrix(0, 0, 6)

 out <- scdd(qux, representation = "H")
 print(out)

 scdd(out$output)

 out <- scdd(qux, representation = "H", adjacency = TRUE,
     inputadjacency = TRUE, incidence = TRUE, inputincidence = TRUE)
 print(out)

 options(error=dump.frames)

 out <- scdd(qux, representation = "V")

