## ------------------------------------------------------------------------
library(spdep)

## ------------------------------------------------------------------------
library(maptools)
columbus <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1])
row.names(columbus)[1:10]

## ------------------------------------------------------------------------
nb_q <- poly2nb(columbus)
nb_q
attr(nb_q, "region.id")[1:10]
is.symmetric.nb(nb_q)

## ------------------------------------------------------------------------
col2 <- droplinks(nb_q, 21)
nb_q[[21]]
col2[[21]]
col2
is.symmetric.nb(col2)
coords <- coordinates(columbus)
plot(nb_q, coords, col="grey")
plot(col2, coords, add=TRUE)

## ------------------------------------------------------------------------
nb_B <- nb2listw(col2, style="B", zero.policy=TRUE)
nb_B$style

## ------------------------------------------------------------------------
B <- as(nb_B, "symmetricMatrix")
all(B == t(B))
str(B)
rownames(B)[1:10]

## ------------------------------------------------------------------------
nb_B1 <- mat2listw(as(B, "dgTMatrix"))
nb_B1$style
all.equal(nb_B1$neighbours, col2, check.attributes=FALSE)
all.equal(attr(nb_B1$neighbours, "region.id"), attr(nb_B$neighbours, "region.id"))

## ------------------------------------------------------------------------
rho <- 0.1
sum(log(1 - rho * eigenw(nb_B)))

## ------------------------------------------------------------------------
n <- nrow(B)
I <- Diagonal(n)
class(I - rho * B)
c(determinant(I - rho * B, logarithm=TRUE)$modulus)

## ------------------------------------------------------------------------
nW <- -B
nChol <- Cholesky(nW, Imult=8)
n * log(rho) + (2 * c(determinant(update(nChol, nW, 1/rho))$modulus))

## ------------------------------------------------------------------------
nb_W <- nb2listw(col2, style="W", zero.policy=TRUE)
W <- as(nb_W, "CsparseMatrix")
str(W)
all(W == t(W))

## ------------------------------------------------------------------------
set.seed(1)
x <- runif(n)
r1 <- as.numeric(W %*% x)
r2 <- lag(nb_W, x, zero.policy=TRUE)
all.equal(r1, r2, check.attributes=FALSE)
plot(x, r1, ylim=c(0,1))
c(x[21], r1[21])

## ------------------------------------------------------------------------
rho <- 0.5
sum(log(1 - rho * eigenw(nb_W)))
class(I - rho * W)
c(determinant(I - rho * W, logarithm=TRUE)$modulus)

## ------------------------------------------------------------------------
LU <- lu(I - rho * W)
sum(log(abs(diag(slot(LU, "U")))))

## ------------------------------------------------------------------------
d <- attr(nb_W$weights, "comp")$d
all.equal(d, card(col2))

## ------------------------------------------------------------------------
dW <- Diagonal(n, d) %*% W
all(dW == t(dW))
isd <- Diagonal(n, 1/sqrt(d))
isd[21,21]
Ws <- as(isd %*% dW %*% isd, "symmetricMatrix")
rowSums(Ws)[21]
class(Ws)
c(determinant(I - rho * Ws, logarithm=TRUE)$modulus)

## ------------------------------------------------------------------------
1/range(eigenw(nb_B))
library(igraph)
f2 <- function(x, extra=NULL) {as.vector(B %*% x)}
1/arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8, which="LA", maxiter=200))$values
1/arpack(f2, sym=TRUE, options=list(n=n, nev=1, ncv=8, which="SA", maxiter=200))$values
#1/arpack(f2, sym=TRUE, options=list(n=n, nev=2, ncv=8, which="BE", maxiter=200))$values
# "BE" gives: At line 558 of file dsaup2.f: Fortran runtime error: 
# Index '9' of dimension 1 of array 'bounds' above upper bound of 8

## ------------------------------------------------------------------------
1/range(eigenw(nb_W))
f2 <- function(x, extra=NULL) {as.vector(W %*% x)}
1/arpack(f2, sym=FALSE, options=list(n=n, nev=1, ncv=8, which="LR", maxiter=200))$values
1/arpack(f2, sym=FALSE, options=list(n=n, nev=1, ncv=8, which="SR", maxiter=200))$values

## ------------------------------------------------------------------------
class(B)
object.size(B)
library(igraph)
g1 <- graph.adjacency(B, mode="undirected")
class(g1)
object.size(g1)

## ------------------------------------------------------------------------
B1 <- get.adjacency(g1)
class(B1)
object.size(B1)
all.equal(B, as(as(B1, "dgTMatrix"), "symmetricMatrix"))

## ------------------------------------------------------------------------
res <- n.comp.nb(col2)
table(res$comp.id)

## ------------------------------------------------------------------------
c1 <- clusters(g1)
c1$no == res$nc
all.equal(c1$membership, res$comp.id)
all.equal(c1$csize, c(table(res$comp.id)), check.attributes=FALSE)

## ------------------------------------------------------------------------
W <- as(nb2listw(col2, style="W", zero.policy=TRUE), "CsparseMatrix")
g1W <- graph.adjacency(W, mode="directed", weighted="W")
c1W <- clusters(g1W)
all.equal(c1W$membership, res$comp.id)

## ------------------------------------------------------------------------
is.connected(g1)
dg1 <- diameter(g1)
dg1
sp_mat <- shortest.paths(g1)
str(sp_mat)

## ------------------------------------------------------------------------
nbl10 <- nblag(col2, maxlag=10)
vals <- sapply(nbl10, function(x) sum(card(x)))
zero <- which(vals == 0)
zero[1]-1

## ------------------------------------------------------------------------
lmat <- lapply(nbl10[1:(zero[1]-1)], nb2mat, style="B", zero.policy=TRUE)
mat <- matrix(0, n, n)
for (i in seq(along=lmat)) mat = mat + i*lmat[[i]]
mat[mat==0] <- Inf
diag(mat) <- 0
all.equal(mat, sp_mat, check.attributes=FALSE)

## ------------------------------------------------------------------------
nb_r <- cell2nb(7, 7, type="rook")
nb_rW <- nb2listw(nb_r, style="W")
spdep:::find_q1_q2(nb_rW)

## ------------------------------------------------------------------------
1/range(Re(eigenw(similar.listw(nb_rW))))

## ------------------------------------------------------------------------
spdep:::find_q1_q2(nb_W)
1/range(Re(eigenw(similar.listw(nb_W))))

