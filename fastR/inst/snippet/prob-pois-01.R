
dpois(0,6/3)                           # 0 customers in 1/3 hour
dpois(2,6/3)                           # 2 customers in 1/3 hour
1-ppois(6,6)                           # more than 6 in 1 hour
dpois(6,6)                             # exactly 6 in 1 hour
ppois(5,6)                             # less than 6 in 1 hour
ppois(30,24) - ppois(19,24)            # 20 to 30 customers in 4 hours
