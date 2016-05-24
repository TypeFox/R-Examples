message("*** withPar() ...")

library("R.devices")

x <- 1
message("x=", x)

withPar({
  layout(1:4)

  x <- c(x, 2)
  message("x=", x)
  stopifnot(all(x == 1:2))

  withPar({
    plot(1:10)
    plot(10:1)
    x <- c(x, 3)
    message("x=", x)
  }, pch=4)
  message("x=", x)
  stopifnot(all(x == 1:3))

  withPar({
    plot(1:10)
    plot(10:1)
    x <- c(x, 4)
    message("x=", x)
  }, pch=0, bg="yellow")
  message("x=", x)
  stopifnot(all(x == 1:4))

  x <- c(x, 5)
  message("x=", x)
  stopifnot(all(x == 1:5))
}, mar=c(2,2,1,1))

message("x=", x)
stopifnot(all(x == 1:5))


# Graphical parameters set "manually" are also reset
opar <- par()
withPar({
  par(pch=4L, lwd=3)
  plot(1:10)
})
stopifnot(all.equal(par(), opar))

message("*** withPar() ... DONE")
