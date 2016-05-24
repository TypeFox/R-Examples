par(mar = c(0.2, 0.2, 0.2, 0.2))
for (n in 41:200) {
    set.seed(711)
    plot.new()
    size = c(replicate(n, 1/rbeta(2, 1.5, 4)))
    center = t(replicate(n, runif(2)))
    center = center[rep(1:n, each = 2), ]
    color = apply(replicate(2 * n, sample(c(0:9, LETTERS[1:6]),
        8, replace = TRUE)), 2, function(x) sprintf("#%s", paste(x,
        collapse = "")))
    points(center, cex = size, pch = rep(20:21, n), col = color)
	box()
    # text(0.5, 0.5, n)
    Sys.sleep(0.4)
}
