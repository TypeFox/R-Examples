library(txtplot)

txtplot(1000000:1000004, 100000*c(1, 2, 3, 2, 1))
txtplot(1000000:1000004, 100000*c(1, 2, 3, 2, 1),
        xlab = "Test 1", ylab = "Test 2")
txtplot(rep(0,10), 0:9)
txtplot(0:9, rep(0,10))
txtplot(rep(10,10), 0:9)
txtplot(0:9, rep(10,10))

txtplot(runif(100), runif(100))
a <- seq(0, 3, length=11)
b <- exp(-a)
txtplot(a, b, pch="o")
txtcurve(x/(x+0.1), 0, 1, xlab = "Dose", ylab="Response")

txtplot(a, b, xlim = c(1,2), ylim = c(0,1))

txtboxplot(rexp(100), xlab = "waiting time")

set.seed(1)
txtplot(runif(10000,0,115000),100001:110000)
txtplot(runif(10000,-115000,0), seq(1,9,length=10000))
txtboxplot(runif(10000,0,115000))
txtboxplot(runif(10000,-115000.0))
