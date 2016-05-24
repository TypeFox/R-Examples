library(RTriangle)
p <- pslg(P=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
          S=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)),
          PA=matrix(1:10, 5, 2))
tp <- triangulate(p, a=0.1)
print(tp)
