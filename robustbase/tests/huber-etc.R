library(robustbase)

### Test sets (all kinds odd/even, constant/regular/outlier)

## n = 0,1,2,3 :
x0 <- numeric(0)
x1 <- 3
x2 <- 1:2
x3 <- c(1:2,10)
## constant (0 mad) + 0--2 outliers
xC <-  rep(1, 12)
xC. <- rep(1, 11)
xC1  <- c(xC,  10)
xC1. <- c(xC., 10)
xC2  <- c(xC1,  100)
xC2. <- c(xC1., 100)
## "uniform"  + 0--2 outliers
y  <- 1:10
y. <- 1:11
y1  <- c(y,  100)
y1. <- c(y., 100)
y2  <- c(y1,  1000)
y2. <- c(y1., 1000)

nms <- ls(pat="^[xy]"); nms; names(nms) <- nms
lx <- lapply(nms,
             function(n) {
                 x <- get(n)
                 m <- mad(x)
                 hx <-
                     if(!is.na(m) && m > 0) MASS::huber(x)
                     else list(m=NA, s=NA)
                 hMx <- huberM(x)
                 list(loc =
                      c(median = median(x),
                        huber  =  hx$m,
                        huberM = hMx$m),
                      scale=
                      c(mad    = m,
                        huber  =  hx$s,
                        huberM = hMx$s))
             })
r <- list(mu = sapply(lx, function(x) x$loc),
          s  = sapply(lx, function(x) x$scale))
r

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
