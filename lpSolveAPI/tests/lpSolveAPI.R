library(lpSolveAPI)

dirSep <- .Platform$file.sep
tests.dir <- paste(system.file(package = "lpSolveAPI"), "tests", sep = dirSep)

################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "absolute1.lp", sep = dirSep))

m <- dim(lp)[1]
n <- dim(lp)[2]

R <- make.lp(0, n)

obj <- numeric(n)
for(j in 1:n) obj[j] <- get.mat(lp, 0, j)
set.objfn(R, obj)

for(i in 1:m) {
  xt <- numeric(n)
  for(j in 1:n) xt[j] <- get.mat(lp, i, j)
  add.constraint(R, xt, get.constr.type(lp)[i], get.constr.value(lp)[i])
}

dimnames(R) <- list(paste("R", 1:9, sep = ""),
                    c("x1", "x2", "x3", "x4", "x2abs", "x4abs"))

lp.control(R, sense = "max")$sense
set.bounds(R, lower = rep(-Inf, 2), upper = rep(Inf, 2), columns = c(2, 4))
set.bounds(R, lower = 1.1, upper = 10, columns = 3)

solve(lp)
solve(R)

vars <- c(3.75, 1.25, 1.1, -0.25, 1.25, 0.25)
opt <- 2.6

stopifnot(all.equal(get.variables(lp), vars))
stopifnot(all.equal(get.variables(R), vars))

stopifnot(all.equal(get.objective(lp), opt))
stopifnot(all.equal(get.objective(R), opt))

rm(lp, m, n, R, obj, j, i, xt, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "binary1.lp", sep = dirSep))

R <- make.lp(4, 4)
dimnames(R) <- list(paste("r", 1:4, sep = "_"), paste("x", 1:4, sep = ""))

set.column(R, 1, c(1, 2, -1, 0))
set.column(R, 2, c(1, -1, 3, 0))
set.column(R, 3, c(0, 0, 0, 1))
set.column(R, 4, c(0, 0, 0, 1))

set.objfn(R, c(-1, -2, 0.1, 3))

set.constr.type(R, c("<=", rep(">=", 3)))
set.rhs(R, c(5, 0, 0, 0.5))

set.type(R, 3, "binary")

solve(lp)
solve(R)

vars <- c(1.66667, 3.33333, 1, 0)
opt <- -8.23333

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "free1.lp", sep = dirSep))
mps <- read.lp(paste(tests.dir, "free1.mps", sep = dirSep))

R <- make.lp(0, 4)
row.add.mode(R, "on")
set.objfn(R, c(1, 2, -4, -3))
add.constraint(R, c(1, 1, 0, 0), "<=", 5)
add.constraint(R, c(2, -1, 0, 0), ">=", 0)
add.constraint(R, c(-1, 3, 0, 0), ">=", 0)
add.constraint(R, c(0, 0, 1, 1), ">=", 0.5)
row.add.mode(R, "off")

dimnames(R) <- list(paste("R", 1:4, sep = ""), paste("x", 1:4, sep = ""))

lp.control(R, sense = "max")$sense
set.bounds(R, lower = rep(-Inf, 2), upper = rep(Inf, 2), columns = c(2, 4))
set.bounds(R, lower = 1.1, upper = 10, columns = 3)

solve(lp)
solve(mps)
solve(R)

vars <- c(1.66667, 3.33333, 1.1, -0.6)
opt <- 5.73333

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(mps), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(-get.objective(mps), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, mps, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "integer1.lp", sep = dirSep))
mps <- read.lp(paste(tests.dir, "integer1.mps", sep = dirSep))

R <- make.lp(4, 4)
dimnames(R) <- list(paste("r", 1:4, sep = "_"), paste("x", 1:4, sep = ""))

set.column(R, 1, c(1, 2, -1, 0))
set.column(R, 2, c(1, -1, 3, 0))
set.column(R, 3, c(0, 0, 0, 1))
set.column(R, 4, c(0, 0, 0, 1))

set.objfn(R, c(-1, -2, 0.1, 3))

set.constr.type(R, c("<=", rep(">=", 3)))
set.rhs(R, c(5, 0, 0, 0.5))

set.type(R, 3, "integer")
set.bounds(R, lower = 1.1, columns = 3)

solve(lp)
solve(mps)
solve(R)

vars <- c(1.66667, 3.33333, 2, 0)
opt <- -8.13333

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(mps), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(mps), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, mps, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "semi-continuous1.lp", sep = dirSep))
mps <- read.lp(paste(tests.dir, "semi-continuous1.mps", sep = dirSep))

R <- make.lp(4, 4)
dimnames(R) <- list(paste("r", 1:4, sep = "_"), paste("x", 1:4, sep = ""))

set.column(R, 1, c(1, 2, -1, 0))
set.column(R, 2, c(1, -1, 3, 0))
set.column(R, 3, c(0, 0, 0, 1))
set.column(R, 4, c(0, 0, 0, 1))

set.objfn(R, c(1, 2, -4, -3))

set.constr.type(R, c("<=", rep(">=", 3)))
set.rhs(R, c(5, 0, 0, 0.5))

lp.control(R, sense = "max")$sense
set.bounds(R, lower = 1.1, upper = 10, columns = 3)
set.semicont(R, 3)

vars <- c(1.66667, 3.33333, 0, 0.5)
opt <- 6.83333

solve(lp)
solve(mps)
solve(R)

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(mps), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(-get.objective(mps), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, mps, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "semi-continuous2.lp", sep = dirSep))

R <- make.lp(4, 4)
dimnames(R) <- list(paste("R", 1:4, sep = ""), paste("x", 1:4, sep = ""))

set.column(R, 1, c(1, 2, -1, 0))
set.column(R, 2, c(1, -1, 3, 0))
set.column(R, 3, c(0, 0, 0, 1))
set.column(R, 4, c(0, 0, 0, 1))

set.objfn(R, c(1, 2, -0.1, -3))

set.constr.type(R, c("<=", rep(">=", 3)))
set.rhs(R, c(5, 0, 0, 0.5))

lp.control(R, sense = "max")$sense
set.bounds(R, lower = 1.1, upper = 10, columns = 3)
set.semicont(R, 3)

vars <- c(1.66667, 3.33333, 1.1, 0)
opt <- 8.22333

solve(lp)
solve(R)

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "semi-continuous3.lp", sep = dirSep))

R <- make.lp(4, 4)
dimnames(R) <- list(paste("R", 1:4, sep = ""), paste("x", 1:4, sep = ""))

set.column(R, 1, c(1, 2, -1, 0))
set.column(R, 2, c(1, -1, 3, 0))
set.column(R, 3, c(0, 0, 0, 1))
set.column(R, 4, c(0, 0, 0, 1))

set.objfn(R, c(1, 2, -0.1, -3))

set.constr.type(R, c("<=", rep(">=", 3)))
set.rhs(R, c(5, 0, 0, 0.5))

lp.control(R, sense = "max")$sense
set.bounds(R, lower = 1.1, upper = 10, columns = 3)
set.semicont(R, 3)
set.type(R, 3, "integer")

vars <- c(1.66667, 3.33333, 2, 0)
opt <- 8.13333

solve(lp)
solve(R)

stopifnot(all.equal(get.variables(lp), vars, tolerance = 1e-5))
stopifnot(all.equal(get.variables(R), vars, tolerance = 1e-5))

stopifnot(all.equal(get.objective(lp), opt, tolerance = 1e-5))
stopifnot(all.equal(get.objective(R), opt, tolerance = 1e-5))

rm(lp, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "sos1.lp", sep = dirSep))
mps <- read.lp(paste(tests.dir, "sos1.mps", sep = dirSep))

R <- make.lp(2, 5)
dimnames(R) <- list(paste("Rc", 1:2, sep = ""), paste("x", 1:5, sep = ""))

set.objfn(R, c(-1, -1, -3, -2, -2))
set.row(R, 1, c(-1, -1, 1, 1, 0))
set.row(R, 2, c(1, 0, 1, -3, 0))
set.constr.type(R, rep("<=", 2))
set.rhs(R, rep(30, 2))

set.bounds(R, upper = c(40, 1, Inf, Inf, 1))
add.SOS(R, "SOS1", 1, 1, 1:5, rep(1, 5))

vars <- c(0, 0, 30, 0, 0)
opt <- -90

solve(lp)
solve(mps)
solve(R)


# These tests fail w/ clang-3.5
#stopifnot(all.equal(get.variables(lp), vars))
#stopifnot(all.equal(get.variables(mps), vars))
#stopifnot(all.equal(get.variables(R), vars))

#stopifnot(all.equal(get.objective(lp), opt))
#stopifnot(all.equal(get.objective(mps), opt))
#stopifnot(all.equal(get.objective(R), opt))

rm(lp, mps, R, vars, opt)


################################################################################
################################################################################

lp <- read.lp(paste(tests.dir, "sos2.lp", sep = dirSep))
mps <- read.lp(paste(tests.dir, "sos2.mps", sep = dirSep))

R <- make.lp(2, 5)
dimnames(R) <- list(paste("Rc", 1:2, sep = ""), paste("x", 1:5, sep = ""))

set.objfn(R, c(-1, -1, -3, -2, -2))
set.row(R, 1, c(-1, -1, 1, 1, 0))
set.row(R, 2, c(1, 0, 1, -3, 0))
set.constr.type(R, rep("<=", 2))
set.rhs(R, rep(30, 2))

set.bounds(R, upper = c(40, 1, Inf, Inf, 1))
add.SOS(R, "SOS2", 2, 1, 1:5, rep(1, 5))

vars <- c(0, 1, 30, 0, 0)
opt <- -91

solve(lp)
solve(mps)
solve(R)

# These tests will also probably fail w/ clang-3.5
#stopifnot(all.equal(get.variables(lp), vars))
#stopifnot(all.equal(get.variables(mps), vars))
#stopifnot(all.equal(get.variables(R), vars))

#stopifnot(all.equal(get.objective(lp), opt))
#stopifnot(all.equal(get.objective(mps), opt))
#stopifnot(all.equal(get.objective(R), opt))

rm(lp, mps, R, vars, opt)


################################################################################
################################################################################
#
#lp <- read.lp(paste(tests.dir, "sos3.lp", sep = dirSep))
#mps <- read.lp(paste(tests.dir, "sos3.mps", sep = dirSep))
#
#R <- make.lp(2, 5)
#dimnames(R) <- list(paste("Rc", 1:2, sep = ""), paste("x", 1:5, sep = ""))
#
#set.objfn(R, c(-1, -1, -3, -2, -2))
#set.row(R, 1, c(-1, -1, 1, 1, 0))
#set.row(R, 2, c(1, 0, 1, -3, 0))
#set.constr.type(R, rep("<=", 2))
#set.rhs(R, rep(30, 2))
#
#set.bounds(R, upper = c(40, 1, Inf, Inf, 1))
#add.SOS(R, "SOS3", 3, 1, 1:5, rep(1, 5))
#
#vars <- c(0, 1, 30.75, 0.25, 0)
#opt <- -93.75
#
#solve(lp)
#solve(mps)
#solve(R)
#
#stopifnot(all.equal(get.variables(lp), vars))
#stopifnot(all.equal(get.variables(mps), vars))
#stopifnot(all.equal(get.variables(R), vars))
#
#stopifnot(all.equal(get.objective(lp), opt))
#stopifnot(all.equal(get.objective(mps), opt))
#stopifnot(all.equal(get.objective(R), opt))
#
#rm(lp, mps, R, vars, opt)
#
#
################################################################################
################################################################################
#
#lp <- read.lp(paste(tests.dir, "sos4.lp", sep = dirSep))
#mps <- read.lp(paste(tests.dir, "sos4.mps", sep = dirSep))
#
#R <- make.lp(2, 5)
#dimnames(R) <- list(paste("Rc", 1:2, sep = ""), paste("x", 1:5, sep = ""))
#
#set.objfn(R, c(-1, -1, -3, -2, -2))
#set.row(R, 1, c(-1, -1, 1, 1, 0))
#set.row(R, 2, c(1, 0, 1, -3, 0))
#set.constr.type(R, rep("<=", 2))
#set.rhs(R, rep(30, 2))
#
#set.bounds(R, upper = c(40, 1, Inf, Inf, 1))
#add.SOS(R, "SOS4", 4, 1, 1:5, rep(1, 5))
#
#vars <- c(40, 1, 50.75, 20.25, 0)
#opt <- -233.75
#
#solve(lp)
#solve(mps)
#solve(R)
#
#stopifnot(all.equal(get.variables(lp), vars))
#stopifnot(all.equal(get.variables(mps), vars))
#stopifnot(all.equal(get.variables(R), vars))
#
#stopifnot(all.equal(get.objective(lp), opt))
#stopifnot(all.equal(get.objective(mps), opt))
#stopifnot(all.equal(get.objective(R), opt))
#
#rm(lp, mps, R, vars, opt)
#
#
################################################################################
################################################################################
#
#lp <- read.lp(paste(tests.dir, "sos5.lp", sep = dirSep))
#mps <- read.lp(paste(tests.dir, "sos5.mps", sep = dirSep))
#
#R <- make.lp(2, 5)
#dimnames(R) <- list(paste("Rc", 1:2, sep = ""), paste("x", 1:5, sep = ""))
#
#set.objfn(R, c(-1, -1, -3, -2, -2))
#set.row(R, 1, c(-1, -1, 1, 1, 0))
#set.row(R, 2, c(1, 0, 1, -3, 0))
#set.constr.type(R, rep("<=", 2))
#set.rhs(R, rep(30, 2))
#
#set.bounds(R, upper = c(40, 1, Inf, Inf, 1))
#add.SOS(R, "SOS5", 5, 1, 1:5, rep(1, 5))
#
#vars <- c(40, 1, 50.75, 20.25, 1)
#opt <- -235.75
#
#solve(lp)
#solve(mps)
#solve(R)
#
#stopifnot(all.equal(get.variables(lp), vars))
#stopifnot(all.equal(get.variables(mps), vars))
#stopifnot(all.equal(get.variables(R), vars))
#
#stopifnot(all.equal(get.objective(lp), opt))
#stopifnot(all.equal(get.objective(mps), opt))
#stopifnot(all.equal(get.objective(R), opt))
#
#rm(lp, mps, R, vars, opt)
#
#
################################################################################
################################################################################

rm(dirSep, tests.dir)


################################################################################
################################################################################

lp <- make.lp(0, 3)

set.objfn(lp, c(1, 1, 1))

add.constraint(lp, c(0, 1, 0), "=", 2)
add.constraint(lp, c(1, 1, 0), "<=", 6)
add.constraint(lp, c(1, 0, -2), ">=", 2)
add.constraint(lp, c(0, 1, 3), "<=", 5)

lp.control(lp, sense = "max")$sense
lp.control(lp, presolve = c("rows"))$presolve

print(lp)
solve(lp)
print(lp)
rm(lp)


################################################################################
################################################################################

lp <- make.lp(0, 3)

set.objfn(lp, c(1, 1, 1))

add.constraint(lp, c(0, 1, 0), "=", 2)
add.constraint(lp, c(1, 1, 0), "<=", 6)
add.constraint(lp, c(1, 0, -2), ">=", 2)
add.constraint(lp, c(0, 1, 3), "<=", 5)

lp.control(lp, sense = "max")$sense
lp.control(lp, presolve = c("cols", "rows"))$presolve

print(lp)
solve(lp)
print(lp)
rm(lp)


################################################################################
################################################################################


