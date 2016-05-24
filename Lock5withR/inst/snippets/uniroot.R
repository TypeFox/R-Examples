# uniroot finds when a function is 0, so we need to build such a function;
# this function is 0 when the margin of error is 0.03:
f <- makeFun(1.96 * sqrt(.5 * .5/ n) - 0.03  ~ n)  
# uniroot needs a function and a lower bound and upper bound to search between
uniroot(f, c(1,50000))$root  

