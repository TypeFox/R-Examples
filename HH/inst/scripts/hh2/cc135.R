
data(cc135)
a1c <-  aov(terms(yield ~ cow + square:period + treat + res.treat,
                  keep.order=TRUE), data=cc135)
summary(a1c)

a1cr <- aov(terms(yield ~ cow + square:period + res.treat + treat,
                  keep.order=TRUE), data=cc135)
summary(a1cr)
