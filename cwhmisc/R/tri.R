"triplot" <- function(a, f, m, symb = 2, grid = FALSE, ...)
{
	ta <- paste(substitute(a))
	tf <- paste(substitute(f))
	tm <- paste(substitute(m))
	tot <- 100/(a + f + m)
	b <- f * tot
	y <- b * .878
	x <- m * tot + b/2
	par(pty = "s")
	oldcol <- par("col")
	plot(x, y, axes = FALSE, xlab = "", ylab = "", xlim = c(-10, 110), ylim
		 = c(-10, 110), type = "n", ...)
        points(x,y,pch=symb)
	par(col = oldcol)
	.trigrid(grid)
	text(-5, -5, ta)
	text(105, -5, tm)
	text(50, 93, tf)
	par(pty = "m")
	invisible()
}

".trigrid" <- function(grid = FALSE)
{
	lines(c(0, 50, 100, 0), c(0, 87.8, 0, 0))	#draw frame
	if(!grid) {
		for(i in 1:4 * 20) {
			lines(c(i, i - 1), c(0, 2 * .878))	#side a-c (base)
			lines(c(i, i + 1), c(0, 2 * .878))
			T.j <- i/2	#side a-b (left)
			lines(c(T.j, T.j + 2), c(i * .878, i * .878))
			lines(c(T.j, T.j + 1), c(i * .878, (i - 2) * .878))
			T.j <- 100 - i/2	#side b-c (right)
			lines(c(T.j, T.j - 2), c(i * .878, i * .878))
			lines(c(T.j, T.j - 1), c(i * .878, (i - 2) * .878))
		}
	}
	else {
		for(i in 1:4 * 20) {
# draw dotted grid
			lines(c(i, i/2), c(0, i * .878), lty = 4, col = 3)	#
			lines(c(i, (50 + i/2)), c(0, .878 * (100 - i)), lty = 4,
				col = 3)	# /
			lines(c(i/2, (100 - i/2)), c(i * .878, i * .878), lty
				 = 4, col = 3)	# -
		}
		par(lty = 1, col = 1)
	}
}
