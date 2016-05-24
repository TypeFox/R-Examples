f1 <- makeFun(qnorm(.99) * sqrt(.5 * .5/ n) - 0.03 ~ n)
uniroot(f1, c(1,50000))$root
f2 <- makeFun(qnorm(.975) * sqrt(.5 * .5/ n) - 0.005 ~ n)
uniroot(f2, c(1,50000))$root
f3 <- makeFun(qnorm(.99) * sqrt(.1 * .9/ n) - 0.005 ~ n)
uniroot(f3, c(1,50000))$root

