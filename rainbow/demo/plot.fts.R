## This demo considers functional data visualization using fts().

library(rainbow)

plot(fts(x = 15:49, y = Australiasmoothfertility$y, xname = "Age", yname = "Fertility rate"))

