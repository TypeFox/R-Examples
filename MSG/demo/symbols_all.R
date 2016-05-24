set.seed(123)
par(mar = c(6, 1, 0.1, 0.1), las = 2)
plot(1:6, rep(1, 6), type = "n", axes = FALSE,
    ann = FALSE, xlim = c(0.8, 7.2), ylim = c(0, 5),
    panel.first = grid())
x = seq(1, 2, length = 4)
symbols(x, 1:4, circles = runif(4, 0.1, 0.6),
    add = TRUE, bg = heat.colors(4), inches = FALSE)
symbols(x + 1, 1:4, squares = runif(4, 0.1,
    0.6), add = TRUE, bg = terrain.colors(4), inches = FALSE)
symbols(x + 2, 1:4, rect = matrix(runif(8,
    0.1, 0.6), 4), add = TRUE, bg = rainbow(4), inches = FALSE)
symbols(x + 3, 1:4, stars = matrix(runif(20,
    0.1, 0.6), 4), add = TRUE, bg = topo.colors(4), inches = FALSE)
symbols(x + 4, 1:4, therm = matrix(runif(12,
    0.1, 0.7), 4), add = TRUE, inches = FALSE)
symbols(x + 5, 1:4, boxplot = matrix(runif(20,
    0.1, 0.7), 4), add = TRUE, inches = FALSE)
axis(1, 1:6, c("circles", "squares", "rectangles",
    "stars", "thermometers", "boxplots"), cex.axis = 0.85)
