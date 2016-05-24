library(ei)
data(sample)
formula = t ~ x
dbuf = ei(formula=formula,total="n",data=sample)
summary(dbuf)
eiread(dbuf, "betab", "betaw")
plot(dbuf, "tomog", "betab", "betaw", "xtfit")
