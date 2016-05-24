# testing beta_1 = 2
t <- ( 1.0671 - 2.0) / 0.0115; t
2 * pt(-abs(t) , df=2361)
# testing beta_1 = 1
t <- ( 1.0671 - 1.0) / 0.0115; t
2 * pt(-abs(t) , df=2361)
# testing beta_2 = 1
t <- ( 0.9206 - 1.0) / 0.0289; t
2 * pt(-abs(t) , df=2361)
