# A. Define the new function of fun()
fun <- function(x, a = 3) {
  y <- (x + a)^2 + 20
  return(y)
}

# B. Supply fun() as an argument to curve(); 1st usage
# win.graph(width = 5.5, height = 2.5, pointsize = 9)
pdf(file = "C:/aErer/fig_optimize.pdf", width = 5.5, height = 2.5,
  pointsize = 9)
par(mai = c(0.4, 0.5, 0.1, 0.1), mgp = c(1.7, 0.7, 0), family = "serif")
curve(expr = fun, from = -15, to = 15, n = 200, ylim = c(-80, 700),
  ylab = expression(y == (x + a)^2 + 20))
curve(expr = fun(x, a = 10), from = -15, to = 15, lty = 2, add = TRUE)
legend(x = "topleft", legend = c("a = 3", "a = 10"),
  lty = c(1, 2), box.lty = 0)  

points(x = c(-10, -3, 15, 15), y = c(20, 20, 344, 645), pch = 19, 
  cex = 1.5, col = "gray40")
text(x = -3, y = -30, labels = "(-3, 20)")
text(x = -10, y = -30, labels = "(-10, 20)")
text(x = 12.5, y = 344, labels = "(15, 344)")
text(x = 12.5, y = 645, labels = "(15, 645)")
dev.off()
    
# C. Find the minimum and maximum values within an interval
# C1. Supply fun() as an argument to optimize(); 2nd usage
ra <- optimize(f = fun, interval = c(-15, 15), maximum = FALSE)
rb <- optimize(f = fun, interval = c(-15, 15), maximum = FALSE, a = 10)
rc <- optimize(f = fun, interval = c(-15, 15), maximum = TRUE)
rd <- optimize(f = fun, interval = c(-15, 15), maximum = TRUE, a = 10)
unlist(ra); unlist(rb); unlist(rc); unlist(rd)

# C2. Supply an anonymous new function as an argument to optimize() 
ra2 <- optimize(f = function(x, a = 3) {(x + a)^2 + 20},
  interval = c(-15, 15), maximum = FALSE)
identical(ra, ra2)  # TRUE