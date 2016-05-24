# Beckett and Diaconis  tack rolling example see Jun Liu (Annals, 1996)
data(tacks)
x <- tacks$x
k <- tacks$k
f <- Bmix(x,k,verb = 5)
plot(f, xlab = "", main = "Beckett-Diaconis Tacky Mixture")
abline(v = c(0.4385, 0.6013, 0.8144),col = 2) # Liu's EM estimates


