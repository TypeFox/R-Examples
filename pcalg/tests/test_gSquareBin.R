library(pcalg)

cat("doExtras:", (doExtras <- pcalg:::doExtras()), "\n")

##################################################
## Simple 2 Variables  x1 -> x2
##################################################
set.seed(123)

## only 2 Variables
b0 <- 0
b1 <- 1
n <- 1000
x1 <- rbinom(n, 1, p = 1/2) ## == sample(c(0,1), n, replace=TRUE)
x2 <- rbinom(n, 1, p <- plogis(b0 + b1* x1))

res1 <- fisher.test(x1,x2)$p.value         # 6.14 e-21
res2 <- gSquareBin(1,2,NULL, cbind(x1,x2)) # 3.74 e-21

(eq1 <- all.equal(c(res1,     res2),
                  c(6.142e-21, 3.7401e-21), tol = 1e-5))
ok1 <- isTRUE(eq1)

##################################################
## chain of 3 variables: x1 -> x2 -> x3
##################################################
b0 <- 0
b1 <- 1
b2 <- 1
n <- 10000

set.seed(12)
x1 <- rbinom(n, 1, p=1/2)
x2 <- rbinom(n, 1, p2 <- plogis(b0 + b1* x1))
x3 <- rbinom(n, 1, p3 <- plogis(b0 + b2* x2))
##
dat <- cbind(x1,x2,x3)
##
system.time(
pv.2 <- c(gSquareBin(3,1, 2, dat),
          gSquareBin(1,3, 2, dat),
          gSquareBin(1,2, 3, dat),
          gSquareBin(3,2, 1, dat))
)

(eq2 <- c(all.equal(pv.2, c(rep(0.89393365, 2), 0, 0), tol=1e-8),
          all.equal(pv.2[1], pv.2[2], tol=1e-15)))

ok2 <- is.logical(eq2) && all(eq2)

##################################################
## collider: x1 -> x3 <- x2
##################################################
b0 <- 0
b1 <- 1.3
b2 <- 1.7
n <- 10000
set.seed(13)
x1 <- rbinom(n, 1, p=1/2)
x2 <- rbinom(n, 1, p=1/2)
x3 <- rbinom(n, 1, p3 <- plogis(b0 + b1*x1 + b2*x2))
##
dat <- cbind(x1,x2,x3)
##
system.time(
pv3 <- c(gSquareBin(1,2, 3,  dat),
         gSquareBin(1,2,NULL,dat),
         gSquareBin(2,1, 3,  dat),
         gSquareBin(1,3, 2,  dat))
)
(eq3 <- all.equal(pv3, c(0, 0.9054922, 0,0), tol=1e-8))
ok3 <- isTRUE(eq3)


##################################################
## add noise variables to collider model
##################################################
x4 <- sample(c(0,1),n, replace=TRUE)
x5 <- sample(c(0,1),n, replace=TRUE)
x6 <- sample(c(0,1),n, replace=TRUE)

dat <- cbind(x1,x2,x3,x4,x5,x6)

system.time(
pv4 <- c(gSquareBin(1,2, c(3,4,5,6), dat),
         gSquareBin(1,3, c(2,4,5,6), dat),
         gSquareBin(2,3, c(1,4,5,6), dat),
         gSquareBin(2,4, c(1,3,5,6), dat),
         gSquareBin(4,2, c(1,3,5,6), dat), # same
         gSquareBin(3,4, c(1,2,5,6), dat),
         ##
         gSquareBin(4,2, c(3,5,6),   dat))
)

(eq4 <- c(all.equal(pv4, c(0, 0, 0, 0.0269101, 0.0269101, 0.7260011, 0.1121826),
                   tol = 2e-7),
          all.equal(pv4[4], pv4[5], tol=1e-15)))

ok4 <- is.logical(eq4) && all(eq4)


##################################################
## rectangle model
##################################################
b0 <- 0
b1 <- 1
b2 <- 1
b3 <- 0.5
b4 <- 0.5
n <- 10000
set.seed(31)
x1 <- rbinom(n, size=1, prob = 1/2)# = sample(c(0,1), n, replace=TRUE)
x2 <- rbinom(n, size=1, prob = (p2 <- plogis(b0 + b1* x1)))
x3 <- rbinom(n, size=1, prob = (p3 <- plogis(b0 + b2* x1)))
x4 <- rbinom(n, size=1, prob = (p4 <- plogis(b0 + b3* x2 + b4* x3)))
##
dat <- cbind(x1,x2,x3,x4)

system.time(
pv5 <- c(gSquareBin(1,4,c(2,3),dat),
         gSquareBin(1,4,  2,   dat),
         gSquareBin(1,4,  3,   dat),
         gSquareBin(1,4,NULL,  dat),
         gSquareBin(2,3,  4,   dat),
         gSquareBin(2,3,  1,   dat),
         gSquareBin(2,3,c(1,4),dat))
)

(eq5 <- all.equal(pv5, c(0.38484, 0.13601, 0.10974, 0, 0.00023, 0.56828, 0.43436), 4e-4))
ok5 <- isTRUE(eq5)

##################################################
## rectangle model -- extended
##################################################
x5 <- rbinom(n, size=1, prob = (p5 <- plogis(b0 + x2)))
x6 <- rbinom(n, size=1, prob = (p6 <- plogis(b0 + x3/2 - x4)))
x7 <- rbinom(n, size=1, prob = (p7 <- plogis(b0 + x5 - x6)))
x8 <- rbinom(n, size=1, prob = (p8 <- plogis(b0 - x3 + x6/2)))
x9 <- rbinom(n, size=1, prob = (p9 <- plogis(b0 + x2 - x7)))
##
dat <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)

## do many tests:
p <- ncol(dat) # = 9
## all (x,y) pairs:
xy <- combn(p,2)
"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
ip <- 1:p
str(rest <- apply(xy, 2, function(e) ip %w/o% e)) # 7 x 36
set.seed(101)
in2 <- seq_len(ncol(xy))
args <- c(lapply(in2, function(i) list(x=xy[1,i], y=xy[2,i], S=rest[,i])),
          lapply(in2, function(i) list(x=xy[1,i], y=xy[2,i], S=sample(rest[,i], size= p-2-1))),
          lapply(in2, function(i) list(x=xy[1,i], y=xy[2,i], S=sample(rest[,i], size= sample(0:(p-2), 1)))),
          lapply(in2, function(i) list(x=xy[1,i], y=xy[2,i], S=sample(rest[,i], size= sample(2:(p-2), 1)))))

if(doExtras) {
    ia <- TRUE
} else {
    ia <- sample(length(args), 25)
    cat("For speed, only testing sub-sample, indices ia =\n"); print(ia)
    args <- args[ia]
}


system.time(
pv6 <- vapply(args, function(L) do.call(gSquareBin, c(L, list(dm=dat))), 1)
)# if(doExtras) about one minute

if(FALSE)
    dput(signif(zapsmall(pv6), 4))
pv6.ok <- c(0, 0, 0.6427, 0.4745, 0.9698, 0.4586, 0.5739, 0.6701, 0.1889,
            0, 0, 0.6476, 0.00601, 0.3586, 0, 0, 0.08916, 0, 0.7642, 0, 0.04021,
            0.4694, 0, 0.8799, 0.3966, 0.358, 0.01514, 0, 0.3615, 0.3244,
            0, 0, 0.6647, 0.6839, 0, 0.4136, 0, 0, 0.7023, 0.303, 0.9189,
            0.5501, 0.7896, 0.8239, 0.04095, 0, 0, 0.5531, 0.6374, 0.755,
            0, 0, 0.629, 0, 0.8377, 0, 0.5749, 0.6346, 0, 0.9513, 0.3138,
            0.6293, 0.000153, 0, 0.7855, 0.1503, 0, 0, 0.2033, 0.9071, 0,
            0.2587, 0, 0, 1.08e-05, 0.0002067, 0.03795, 0.5501, 3.5e-06,
            0.001827, 0.5683, 0, 0, 0.005468, 0.4085, 0.1476, 0, 0, 0.6822,
            0, 0.7642, 0, 0.7184, 0.4596, 0, 0.8799, 0.02051, 0.03413, 0.01689,
            0, 0.2786, 0.07084, 0, 0, 5.17e-05, 0.4915, 0, 0.5029, 0, 0,
            0.1891, 0.09992, 0.9698, 0.3937, 5.94e-05, 0.5473, 0.3894, 0,
            0, 0.2236, 2.3e-06, 0.784, 0, 0, 0.3504, 0, 0.7642, 0, 0.3572,
            0.4259, 0, 0.8139, 0.3326, 0.01518, 0.01514, 0, 0.3541, 0.021,
            0, 0, 0.6647, 0.6167, 0, 0.4136)

(eq6 <- all.equal(pv6, pv6.ok[ia], tol= 1e-4))
ok6 <- isTRUE(eq6)

if (!all(ok1,ok2,ok3,ok4,ok5,ok6))
    stop("Test gSquareBin wrong: Some dependence was not estimated correctly!")
