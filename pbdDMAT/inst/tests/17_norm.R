# ################################################
# ------------------------------------------------
# Norm stuff
# ------------------------------------------------
# ################################################


suppressPackageStartupMessages(library(pbdDMAT, quiet=T))

init.grid()

BL <- 4

mn <- 100
sd <- 10

seed <- allreduce(sample(1:5000, size=1))
comm.set.seed(seed=seed, diff=F)

tol <- 1e-8


f1 <- function(type)
{
  out1 <- norm(x, type=type)
  out2 <- norm(x=dx, type=type)
  comm.print(all.equal(out1, out2), quiet=T)
}

f2 <- function(type)
{
    out1 <- rcond(x, norm=type)
    out2 <- rcond(x=dx, norm=type)
    comm.print(all.equal(out1, out2), quiet=T)
}

### Don't check these because R's kappa() is wrong as all hell
#f3 <- function(type)
#{
#    out1 <- kappa(x, norm=type)
#    out2 <- kappa(dx, norm=type)
#    comm.print(all.equal(out1, out2), quiet=T)
#}

# ---------------------------------------------------
# Tests
# ---------------------------------------------------

#######
# m=n #
#######
comm.cat("-----------M=N-----------\n", quiet=T)
M <- 1e2
N <- 1e2
x <- matrix(rnorm(M*N, mean=mn, sd=sd), M, N)
dx <- as.ddmatrix(x, BL)

comm.print("Norms", quiet=T)
f1("M")
f1("F")
f1("O")
f1("I")

comm.print("rcond", quiet=T)
f2("1")
f2("O")
f2("I")

#comm.print("kappa", quiet=T)
#f3("O")
#f3("I")

#####
# m>n #
#######
comm.cat("-----------M>N-----------\n", quiet=T)
M <- 1e3
N <- 1e2
x <- matrix(rnorm(M*N, mean=mn, sd=sd), M, N)
dx <- as.ddmatrix(x, BL)

comm.print("Norms", quiet=T)
f1("M")
f1("F")
f1("O")
f1("I")

comm.print("rcond", quiet=T)
f2("1")
f2("O")
f2("I")

#comm.print("kappa", quiet=T)
#f3("O")
#f3("I")

#######
# m<n #
#######
comm.cat("-----------M<N-----------\n", quiet=T)
M <- 1e2
N <- 1e3
x <- matrix(rnorm(M*N, mean=mn, sd=sd), M, N)
dx <- as.ddmatrix(x, BL)

comm.print("Norms", quiet=T)
f1("M")
f1("F")
f1("O")
f1("I")

comm.print("rcond", quiet=T)
f2("1")
f2("O")
f2("I")

#comm.print("kappa", quiet=T)
#f3("O")
#f3("I")



finalize()
