# earth.times.R

library(earth)
data(ozone1)
library(mda) # for mars

test <- function(nk, degree, use.ozone) {
    if (use.ozone) {
        # example of a small model: 9 predictors 330 cases

        x = ozone1[,-1]
        y = ozone1[,1]
        niters <- 100
    } else {
        # example of a bigger model: 30 predictors 3000 cases

        robot.arm <- function(x) { # lifted from Friedman's Fast MARS paper
          x. <- with(x, l1 * cos(theta1) - l2 * cos(theta1 + theta2) * cos(phi))
          y  <- with(x, l1 * sin(theta1) - l2 * sin(theta1 + theta2) * cos(phi))
          z  <- with(x, l2 * sin(theta2) * sin(phi))
          sqrt(x.^2 + y^2 + z^2)
        }
        set.seed(1)   # for reproducibility
        ncases <- 3000
        l1     <- runif(ncases, 0, 1)
        l2     <- runif(ncases, 0, 1)
        theta1 <- runif(ncases, 0, 2 * pi)
        theta2 <- runif(ncases, 0, 2 * pi)
        phi    <- runif(ncases, -pi/2, pi/2)
        x <- cbind(l1, l2, theta1, theta2, phi)
        for (i in 1:25) # 25 dummy vars, so 30 vars in total
           x <- cbind(x, runif(ncases, 0, 1))
        x <- data.frame(x)
        y <- robot.arm(x)
        niters <- 5
    }
    e.time <- system.time(e <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=0))
    gcv.null <- e$gcv.per.subset[1]
    e.grsq <- 1 - e$gcv/gcv.null

    e.minspan1.time <- system.time(e.minspan1 <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=1))

    e.fastk0.time <- system.time(e.fastk0 <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=0, fast.k=0))

    e.fastk10.time <- system.time(e.fastk10 <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=0, fast.k=10))

    e.fastk5.time <- system.time(e.fastk5 <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=0, fast.k=5))

    e.nobcache.time <- system.time(e.nobcache <- for (i in 1:niters)
        earth(x, y, nk=nk, degree=degree, minspan=0, Use.beta.cache=FALSE))

    # dummy func to estimate time taken by an "allowed" function
    allowed <- function(degree, pred, parents)
    {
        if (degree > 0 && (parents[1] == 999 || pred == 999))
            return(FALSE) # never get here
        TRUE
    }
    e.allowed.time <- system.time(e.allowed <- for (i in 1:niters)
      earth(x, y, nk=nk, degree=degree, minspan=0, allowed=allowed))

    m.time <- system.time(m <- for (i in 1:niters)
         mars(x, y, degree=degree,  nk=nk))
    m.grsq <- 1 - m$gcv/gcv.null

    options(digits=2)
    cat("\t    ",
        nk, "\t",
        degree, "\t",
        e.time[1] / niters, "\t",
        m.time[1]          / e.time[1], "\t",
        e.minspan1.time[1] / e.time[1], "\t\t",
        e.fastk0.time[1]   / e.time[1], "\t",
        e.fastk10.time[1]  / e.time[1], "\t",
        e.fastk5.time[1]   / e.time[1], "\t",
        e.nobcache.time[1] / e.time[1], "\t",
        e.allowed.time[1]  / e.time[1], "\t ",
        e.grsq, "\t",
        m.grsq, "\t",
        e.minspan1$grsq, "\t",
        e.fastk0$grsq, "\t",
        e.fastk10$grsq, "\t",
        e.fastk5$grsq, "\n")
}
cat("data         nk     degree       earth   ")
cat("------------------time ratios---------------------------------  ")
cat("------------------grsq--------------------------\n")
cat("                                 time    ")
cat("mars   minspan=1   no-fastk fastk=10 fastk=5 no-bcache allowed  ")
cat("earth   mars minspan=1 no-fastk fastk=10 fastk=5\n")

cat("ozone 9x330\n");
test(nk=21, degree=1,  use.ozone=1)
test(nk=21, degree=2,  use.ozone=1)
test(nk=21, degree=3,  use.ozone=1)
test(nk=21, degree=10, use.ozone=1)
cat("\n")

test(nk=51, degree=1,  use.ozone=1)
test(nk=51, degree=2,  use.ozone=1)
test(nk=51, degree=3,  use.ozone=1)
test(nk=51, degree=10, use.ozone=1)
cat("\n")

test(nk=101, degree=1,  use.ozone=1)
test(nk=101, degree=2,  use.ozone=1)
test(nk=101, degree=3,  use.ozone=1)
test(nk=101, degree=10, use.ozone=1)

cat("\nrobot arm 30x3000\n");
test(nk=21, degree=1,  use.ozone=0)
test(nk=21, degree=2,  use.ozone=0)
test(nk=21, degree=3,  use.ozone=0)
test(nk=21, degree=10, use.ozone=0)
cat("\n")

test(nk=51, degree=1,  use.ozone=0)
test(nk=51, degree=2,  use.ozone=0)
test(nk=51, degree=3,  use.ozone=0)
test(nk=51, degree=10, use.ozone=0)
cat("\n")

test(nk=101, degree=1,  use.ozone=0)
test(nk=101, degree=2,  use.ozone=0)
test(nk=101, degree=3,  use.ozone=0)
test(nk=101, degree=10, use.ozone=0)
cat("\n")

test(nk=201, degree=1,  use.ozone=0)
test(nk=201, degree=2,  use.ozone=0)
test(nk=201, degree=3,  use.ozone=0)
test(nk=201, degree=10, use.ozone=0)
