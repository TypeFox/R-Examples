1 - (1-p)^(1/6) -> q; q         # P(damage to particular O-ring)
1 - dbinom(0,6,q)               # P(damage to >0 O-rings)
cbind(0:6, dbinom(0:6,6,q))     # table of all probabilities
