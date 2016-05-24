library("pdc")

stopifnot(entropyHeuristic(X = sin(1:1000)^2)$m==5)
