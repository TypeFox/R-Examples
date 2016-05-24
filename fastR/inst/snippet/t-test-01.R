t <- ( 10.3 - 10 )/ (0.4 / sqrt(10)); t;       # test statistic
2 * pt(-abs(t),df=9);            # p-value using t-distribution
2 * pnorm(-abs(t));              # "p-value" using normal distribution
