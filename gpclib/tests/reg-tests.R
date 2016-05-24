suppressMessages(library(gpclib))

## Make some random polygons
set.seed(100)
a <- cbind(rnorm(100), rnorm(100))
a <- a[chull(a), ]

## Convert `a' from matrix to "gpc.poly"
a <- as(a, "gpc.poly")
show(a)

b <- cbind(rnorm(100), rnorm(100))
b <- as(b[chull(b), ], "gpc.poly")
show(b)

## More complex polygons with an intersection
p1 <- read.polyfile(system.file("poly-ex/ex-poly1.txt", package = "gpclib"))
p2 <- read.polyfile(system.file("poly-ex/ex-poly2.txt", package = "gpclib"))

## Plot both polygons and highlight their intersection in red
plot(app <- append.poly(p1, p2))
show(app)
plot(int <- intersect(p1, p2), poly.args = list(col = 2), add = TRUE)
show(int)

## Highlight the difference p1 \ p2 in green
plot(sdif <- setdiff(p1, p2), poly.args = list(col = 3), add = TRUE)
show(sdif)

## Highlight the difference p2 \ p1 in blue
plot(sdif <- setdiff(p2, p1), poly.args = list(col = 4), add = TRUE)
str(sdif)

## Plot the union of the two polygons
plot(un <- union(p1, p2))
str(un)

## Take the non-intersect portions and create a new polygon
## combining the two contours
p.comb <- append.poly(setdiff(p1, p2), setdiff(p2, p1))
str(p.comb)



## Coerce from a matrix
x <- 
structure(c(0.0934073560027759, 0.192713393476752, 0.410062456627342, 
0.470020818875781, 0.41380985426787, 0.271408743927828, 0.100902151283831, 
0.0465648854961832, 0.63981588032221, 0.772382048331416,
0.753739930955121, 0.637744533947066, 0.455466052934407,
0.335327963176065, 0.399539700805524, 
0.600460299194476), .Dim = c(8, 2))
y <- 
structure(c(0.404441360166551, 0.338861901457321, 0.301387925052047, 
0.404441360166551, 0.531852879944483, 0.60117973629424, 0.625537820957668, 
0.179976985040276, 0.341542002301496, 0.445109321058688,
0.610817031070196, 0.596317606444189, 0.459608745684695,
0.215189873417722), .Dim = c(7, 2))

x1 <- as(x, "gpc.poly")
y1 <- as(y, "gpc.poly")

plot(append.poly(x1, y1))
plot(intersect(x1, y1), poly.args = list(col = 2), add = TRUE)


## Show the triangulation
plot(append.poly(x1, y1))
triangles <- triangulate(append.poly(x1,y1))
print(triangles)
