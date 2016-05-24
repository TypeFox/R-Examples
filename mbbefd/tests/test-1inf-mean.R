require(mbbefd)

#check of expectation for oi distribution
n <- 1e4

probs <- c(1/(2:9))

sapply(probs, function(p) 
{
  x <- roiunif(n, p)
c(mean(x), mbbefd:::tmean1(doiunif, p1=p),
  mbbefd:::tmean2(poiunif, p1=p), mbbefd:::tmean3(ecoiunif, p1=p))
}
)



sapply(probs, function(p) 
{
  x <- roistpareto(n, a=2, p)
  c(mean(x), mbbefd:::tmean1(doistpareto, a=2, p1=p),
    mbbefd:::tmean2(poistpareto, a=2, p1=p), mbbefd:::tmean3(ecoistpareto, a=2, p1=p))
}
)




sapply(probs, function(p) 
{
  x <- roibeta(n, shape1=2, shape2=3, p)
  c(mean(x), mbbefd:::tmean1(doibeta, shape1=2, shape2=3, p1=p),
    mbbefd:::tmean2(poibeta, shape1=2, shape2=3, p1=p), 
    mbbefd:::tmean3(ecoibeta, shape1=2, shape2=3, p1=p))
}
)


sapply(probs, function(p) 
{
  x <- roigbeta(n, shape0=pi, shape1=2, shape2=3, p)
  c(mean(x), mbbefd:::tmean1(doigbeta, shape0=pi, shape1=2, shape2=3, p1=p),
    mbbefd:::tmean2(poigbeta, shape0=pi, shape1=2, shape2=3, p1=p), 
    mbbefd:::tmean3(ecoigbeta, shape0=pi, shape1=2, shape2=3, p1=p))
}
)


  x <- rmbbefd(n, a=2, b=1/2)
  c(mean(x), mbbefd:::tmean1(dmbbefd, a=2, b=1/2),
    mbbefd:::tmean2(pmbbefd, a=2, b=1/2), 
    mbbefd:::tmean3(ecmbbefd, a=2, b=1/2))


x <- rmbbefd(n, a=-1/2, b=3)
c(mean(x), mbbefd:::tmean1(dmbbefd, a=-1/2, b=3),
  mbbefd:::tmean2(pmbbefd, a=-1/2, b=3), 
  mbbefd:::tmean3(ecmbbefd, a=-1/2, b=3))
