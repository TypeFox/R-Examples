library(grid);
library(mvtnorm);
set.seed(123);
sigma <- matrix(c(5,2,2,3), ncol = 2);
x <- rmvnorm(n = 200, mean = c(1, 2), sigma = sigma);
range.v1 <- range(x[, 1]);
range.v2 <- range(x[, 2]);
ext.r <- max(diff(range.v1), diff(range.v2)) / 2;
scale.v1 <- mean(range.v1) + c(-1.05, 1.05) * ext.r;
scale.v2 <- mean(range.v2) + c(-1.05, 1.05) * ext.r;

vp <- viewport(x = 0.58, y = 0.58, width = 0.8,
	     	   height = 0.8,
               xscale = scale.v1,
			   yscale = scale.v2);
grid.newpage();
pushViewport(vp);
grid.rect();
grid.xaxis();
grid.yaxis();
grid.text("V1", y = unit(-3, "lines"));
grid.text("V2", x = unit(-3, "lines"), rot = 90);
grid.points(x[, 1], x[, 2], pch = 19,
			gp = gpar(col = "blue", alpha = 0.3));
grid.points(mean(x[, 1]), mean(x[, 2]), pch = 19,
			gp = gpar(col = "red"));
grid.clip();

pr <- princomp(x);
ran.pr <- apply(pr$scores, 2, function(x) max(abs(x)));
scale.pr <- rbind(-1.5 * ran.pr, 1.5 * ran.pr);
size.pr <- 3 * ran.pr / c(diff(scale.v1), diff(scale.v2));
angle.pr <- asin(pr$loadings[1, 2]) / pi * 180;
vp.xaxis <- viewport(x = unit(mean(x[, 1]), "native"),
					 y = unit(mean(x[, 2]), "native"),
					 width = size.pr[1],
					 just = c("center", "bottom"),
					 angle = angle.pr,
					 xscale = scale.pr[, 1],
					 gp = gpar(col = "red"));
pushViewport(vp.xaxis);
grid.lines(c(-2, 2), c(0, 0));
grid.xaxis();
popViewport();
vp.yaxis <- viewport(x = unit(mean(x[, 1]), "native"),
					 y = unit(mean(x[, 2]), "native"),
					 height = size.pr[2],
					 just = c("left", "center"),
					 angle = angle.pr,
					 yscale = scale.pr[, 2],
					 gp = gpar(col = "red"));
pushViewport(vp.yaxis);
grid.lines(c(0, 0), c(-2, 2));
grid.yaxis();

