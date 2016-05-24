m <- 5/(2+5); m                                       # mean
v <- 5 * 2 / ( (5+2)^2 * (5 + 2 + 1) ); v             # var
qbeta(0.5,5,2)                                        # median
pbeta(m,5,2)                                          # less than mean
pbeta(0.4,5,2) - pbeta(0.2,5,2)                       # in a range
pbeta(m+sqrt(v), 5, 2) - pbeta(m-sqrt(v), 5, 2)       # near mean
