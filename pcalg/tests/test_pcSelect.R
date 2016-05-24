#### Testing  pcSelect()  _and_  (its main helper) mcor()
####          ~~~~~~~~                             ~~~~

library(pcalg)

p <- 10
n <- 1000

set.seed(101)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
res1 <- pcSelect(y,dm, alpha=0.05)
if (!all(res1$G == 1:9 %in% c(4,5,6)))
  stop("Test of pcSelect: Consistency problem 101")
dput(signif(res1$zMin, 7))
zM.1 <- c(0.07643311, 0.6568222, 0.7363025,
          11.62001, 6.621493, 18.64382,
          1.447, 1.212824, 1.221812)
stopifnot(all.equal(res1$zMin, zM.1, tol = 1e-6))

## Now all other methods:

(corMeths <- eval(formals(mcor)$method))
## the other methods (but the default):
cMeths <- corMeths[-1]
C. <- zMin. <- setNames(as.list(cMeths), cMeths)
Cstats <- function(C) {
    ## numeric symmetric matrix
    stopifnot(is.matrix(C), is.numeric(C), isSymmetric(C))
    cbind(diag = diag(C),
          colSums = colSums(C),
          e.values = eigen(C, only.values=TRUE)$values)
}

for(meth in cMeths) {
    cat(meth,":\n")
    rr <- pcSelect(y,dm, alpha=0.05, corMethod = meth)
    iG <- which(as.vector(rr$G))
    if (!identical(iG, c(4L,5L,6L)))
        stop(sprintf("pcSelect(dm101.., corMethod=\"%s\"): which(G) = %s\n",
                     meth, paste(iG, collapse=", ")))
    zMin.[[meth]] <- rr$zMin
    C.   [[meth]] <- mcor(d.mat, method=meth)
}

C.st <- lapply(C., Cstats)
stopifnot( sapply(C.st, function(st) st[,"diag"] == 1) )

C.EVal <- sapply(C.st, function(st) st[,"e.values"])
C.sums <- sapply(C.st, function(st) st[,"colSums"])
dns <- list(NULL, c("Qn", "QnStable", "ogkScaleTau2", "ogkQn", "shrink"))

## dput(lapply(zMin., signif, digits=7))
zMin1.Exp <- list(
    Qn = c(0.5316588, 1.088152, 0.8118473, 11.34048,
           6.416001, 18.79211, 0.9179979, 0.4779155, 1.346989),
    QnStable = c(0.9129029, 1.710694, 0.8471604, 10.56123, 6.168661,
                 18.88213, 0.9710747, 0.4337806, 1.140709),
    ogkScaleTau2 = c(0.4698467, 0.886344, 0.5731621, 10.94707, 6.335696,
                     18.60215, 0.9969964, 1.045772, 1.244915),
    ogkQn = c(0.7588748, 0.6270149, 0.6489765, 11.19786, 6.052107,
              18.683, 0.8969601, 0.9808662, 1.118423),
    shrink = c(0.07570146, 0.6493067, 0.7303193, 11.49862, 6.556261,
               18.45523, 1.433129, 1.199766, 1.211877))
##
dput(signif(C.EVal, digits=7))
C.EV1.Exp <- structure(
    c(2.902214, 1.991134, 1.352243, 0.9937615, 0.6887558,
      0.675548, 0.5645716, 0.3164719, 0.3107822, 0.2045182,
      2.895631, 1.980508, 1.3593, 0.9949763, 0.6918006,
      0.6710074, 0.5672133, 0.3192135, 0.3023458, 0.2180048,
      2.885849, 1.966942, 1.340686, 1.003727, 0.7069337,
      0.6863926, 0.5652152, 0.3226591, 0.3090602, 0.2125357,
      2.863277, 1.9834, 1.360826, 0.9985534, 0.70277,
      0.69042, 0.5662882, 0.3179164, 0.3060326, 0.2105163,
      2.884824, 1.989006, 1.336554, 1.000628, 0.7192748,
      0.6763509, 0.5504135, 0.3225079, 0.3041295, 0.216311),
    .Dim = c(10L, 5L), .Dimnames = dns)
##
dput(signif(C.sums, digits=7))
C.sum1.Exp <- structure(
    c(2.717515, 3.48909, 2.081151, 3.301211, 2.762732,
      1.388646, 2.721029, 1.384193, 1.71164, 3.184033, 2.712734, 3.508554,
      2.086835, 3.274161, 2.775645, 1.390581, 2.705537, 1.392348, 1.705256,
      3.164881, 2.744097, 3.500246, 2.132502, 3.314979, 2.721613, 1.407829,
      2.675133, 1.376944, 1.696778, 3.147149, 2.720277, 3.471509, 2.125648,
      3.28348, 2.673719, 1.392961, 2.662582, 1.360947, 1.706638, 3.115694,
      2.790109, 3.431366, 2.075317, 3.269821, 2.714512, 1.462898, 2.714256,
      1.337588, 1.736093, 3.192173),
    .Dim = c(10L, 5L), .Dimnames = dns)

stopifnot(all.equal(zMin. , zMin1.Exp, tol = 1e-6),
          all.equal(C.EVal, C.EV1.Exp , tol = 1e-6),
          all.equal(C.sums, C.sum1.Exp, tol = 1e-6),
          TRUE)

###--------------------------------------- Ex. 2 -------------------------------------

set.seed(456)
myDAG <- randomDAG(p, prob = 0.2)
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
y <- d.mat[,10]
dm <- d.mat[,-10]
iG.Exp <- c(5L,8L,9L)
res2 <- pcSelect(y,dm,alpha=0.05)
if (!identical(which(as.vector(res2$G)), iG.Exp))
  stop("Test of pcSelect: Consistency problem 456")
## dput(signif(res2$zMin, 7))
zM.2 <- c(0.1469954, 0.915466, 0.5442373,
          0.1089084, 14.16709, 1.056202,
          0.4840786, 17.10435, 14.51551)
stopifnot(all.equal(res2$zMin, zM.2, tol = 1e-6))

## Now all other methods:
for(meth in cMeths) {
    cat(meth,":\n")
    rr <- pcSelect(y,dm, alpha=0.05, corMethod = meth)
    iG <- which(as.vector(rr$G))
    if (!identical(iG, iG.Exp))
        (if(all(iG.Exp %in% iG)) message else stop)(
            sprintf("pcSelect(dm456.., corMethod=\"%s\"): which(G) = %s\n",
                    meth, paste(iG, collapse=", ")))
    zMin.[[meth]] <- rr$zMin
    C.   [[meth]] <- mcor(d.mat, method=meth)
}

C.st <- lapply(C., Cstats)
stopifnot( sapply(C.st, function(st) st[,"diag"] == 1) )

C.EVal <- sapply(C.st, function(st) st[,"e.values"])
C.sums <- sapply(C.st, function(st) st[,"colSums"])

##
dput(signif(C.EVal, digits=7))
C.EV2.Exp <- structure(
    c(2.395184, 1.686566, 1.280328, 1.094699, 1.002715,
      0.8997931,0.7362697,0.5018258,0.2396817,0.1629375,
      2.414695, 1.669223, 1.277225, 1.096991, 1.00195,
      0.8970272,0.7380218,0.504836, 0.2350952,0.1649354,
      2.386837, 1.679313, 1.309865, 1.101403, 1.002767,
      0.865312, 0.715178, 0.5221997,0.2352581,0.1818661,
      2.370872, 1.681723, 1.315575, 1.112341, 1.008541,
      0.8657987,0.7116331,0.521083, 0.2312147,0.1812185,
      2.413167, 1.679511, 1.30454,  1.090105, 1.000066,
      0.879334, 0.6960254,0.5191567,0.2335488,0.1845462),
    .Dim = c(10L, 5L), .Dimnames = dns)
##
dput(signif(C.sums, digits=7))
C.sum2.Exp <- structure(
    c(0.9910661,1.77166,  1.919368, 1.752008, 2.239346,
      2.056446, 1.496665, 2.717977, 1.909636, 3.079426,
      0.9969023,1.797362, 1.928749, 1.756807, 2.228714,
      2.072544, 1.517333, 2.795068, 1.921609, 3.127389,
      0.9945908,1.709868, 1.945329, 1.744885, 2.179677,
      2.035526, 1.5155,   2.780851, 1.868081, 3.08998,
      1.00944,  1.692255, 1.943654, 1.754006, 2.144725,
      2.019437, 1.490262, 2.764678, 1.846422, 3.078384,
      1.023379, 1.778081, 2.029191, 1.678237, 2.140316,
      2.036679, 1.470303, 2.754962, 1.936074, 3.135426),
    .Dim = c(10L, 5L), .Dimnames = dns)


## dput(lapply(zMin., signif, digits=7))
zMin2.Exp <-
    list(
        Qn = c(1.055523, 2.555242, 1.691455, 0.7796301, 11.39895,
               0.8692785, 0.04107331, 14.0156, 14.82198),
        QnStable = c(0.6208816, 2.475236, 1.831045, 0.9328639, 11.11562,
                     1.139745, 0.04850462, 14.21195, 14.77496),
        ogkScaleTau2 = c(0.4164263, 1.092748, 1.689739, 0.4691277, 11.92243,
                         1.028907, 0.06462646, 13.80103, 14.60776),
        ogkQn = c(0.237811, 1.916074, 1.704892, 0.3637768, 14.42101,
                  0.8777305, 0.05152724, 16.72384, 14.60674),
        shrink = c(0.1459792, 0.902588, 0.540349, 0.1081555, 13.98006,
                   1.047259, 0.4807316, 16.88456, 14.4006))

stopifnot(all.equal(zMin. , zMin2.Exp , tol = 1e-6),
          all.equal(C.EVal, C.EV2.Exp , tol = 1e-6),
          all.equal(C.sums, C.sum2.Exp, tol = 1e-6),
          TRUE)
