O
a <- sum(O[1:2]) / (2 * sum(O)); a
b <- sum(O[3:4]) / (2 * sum(O)); b
a+b                                  # should equal 0.5
lnum <- 275 * log(a) + 211 * log(b)
ldenom <- sum( O * log (O/ sum(O)))
G <- -2 * ( lnum - ldenom); G
1 - pchisq(G,df=2)
