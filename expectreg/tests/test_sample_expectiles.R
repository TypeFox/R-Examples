require("expectreg")

set.seed(9484470)

#expectile
r <-rnorm(100000)
e <- expectile(r,seq(0.01, 0.99, 0.1))
en <- enorm(seq(0.01, 0.99, 0.1))
stopifnot(max(abs(e-en)) < 0.1)


#pemq, demq, remq, qemq, eemq

y <- seq(.01,.99,.01)
stopifnot(max(abs(eemq(y)-qemq(y)))<.000001)

e <- eemq(y)
w <- pemq(e)
stopifnot(0<w && w<1)
stopifnot(max(abs(y-w))<.0000001)

r <- remq(1000)
w <- pemq(r)
q <- qemq(w)
stopifnot(max(abs(r-q))<.0000001)


#enorm,penorm,ebeta, pebeta, eunif, peunif, et, pet, elnorm, pelnorm, egamma,....

stopifnot(max(abs(penorm(enorm(y))-y))<.000001)
stopifnot(0<penorm(enorm(y))&& penorm(enorm(y))<1)

stopifnot(max(abs(pebeta(ebeta(y))-y))<.000001)
stopifnot(0<pebeta(ebeta(y))&& pebeta(ebeta(y))<1)

stopifnot(max(abs(peunif(eunif(y))-y))<.000001)
stopifnot(0<peunif(eunif(y))&& peunif(eunif(y))<1)

stopifnot(max(abs(pet(et(y,2),2)-y))<.000001)
stopifnot(0<pet(et(y,2),2)&& pet(et(y,2),2)<1)

stopifnot(max(abs(pelnorm(elnorm(y))-y))<.000001)
stopifnot(0<pelnorm(elnorm(y))&& pelnorm(elnorm(y))<1)

stopifnot(max(abs(pegamma(egamma(y,2),2)-y))<.000001)
stopifnot(0<pegamma(egamma(y,2),2)&& pegamma(egamma(y,2),2)<1)

stopifnot(max(abs(peexp(eexp(y))-y))<.000001)
stopifnot(0<peexp(eexp(y))&& peexp(eexp(y))<1)

stopifnot(max(abs(pechisq(echisq(y,2),2)-y))<.000001)
stopifnot(0<pechisq(echisq(y,2),2)&& pechisq(echisq(y,2),2)<1)






