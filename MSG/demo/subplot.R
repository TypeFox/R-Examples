library(grid);
library(KernSmooth);
set.seed(123);
n <- 200;
x <- runif(n);
y <- sin(8 * x) / x + rnorm(n);
scale.x <- extendrange(range(x));
scale.y <- extendrange(range(y));
x0 <- seq(min(x), max(x), length.out = 200);
y0t <- sin(8 * x0) / x0;
y0f <- locpoly(x, y, bandwidth = 0.15, gridsize = 200,
			   range.x = range(x))$y;
vp <- viewport(x = 0.52, y = 0.60, width = 0.8,
	     	   height = 0.77,
               xscale = scale.x,
			   yscale = scale.y);
grid.newpage();
pushViewport(vp);
doPlot <- function(x, y, x0, y0t, y0f)
{
	grid.points(x, y, gp = gpar(cex = 0.5));
	grid.lines(x0, y0t, default.units = "native",
			   gp = gpar(col = "blue"));
	grid.lines(x0, y0f, default.units = "native",
			   gp = gpar(col = "red"));
}
grid.rect();
grid.xaxis();
grid.yaxis();
grid.text("x", y = unit(-3, "lines"));
grid.text("y", x = unit(-3, "lines"), rot = 90);
grid.rect(0.6, -1, 0.4, 2, default.units = "native",
		  gp = gpar(col = NA, fill = rgb(1, 1, 0, 0.8)));
doPlot(x, y, x0, y0t, y0f);

vp.sub <- viewport(x = 0.68, y = 7, width = 0.4 * 1.5,
				   height = 4 * 1.5,
				   default.units = "native",
				   xscale = c(0.4, 0.8),
				   yscale = c(-2, 0));
pushViewport(vp.sub);
grid.rect(gp = gpar(lty = 2));
grid.clip();
doPlot(x, y, x0, y0t, y0f);

