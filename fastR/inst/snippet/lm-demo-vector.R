v0 <- rep(1,4); v0
u0 <- v0/vlength(v0); u0
v1 <- x - mean(x); v1
u1 <- v1/vlength(v1); u1
#
# projecting into the model space
project(y,v0)
project(y,v1)
#
# two ways to compute beta_1-hat
project(y,v1,type='l')/vlength(v1)
dot(y,v1)/(vlength(v1))^2 -> b1; b1
#
# two ways to compute alpha_0-hat
project(y,v0,type='l')/vlength(v0)
dot(y,v0)/(vlength(v0))^2 -> a0; a0
#
# beta_1-hat
a0 - b1 * mean(x)
