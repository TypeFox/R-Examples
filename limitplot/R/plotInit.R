#
# Limitplot functions -- plotInit
#

plotInit <- function(pl, lod, vnum, stack, ratio, axis, main, xlab, ylab, jitterwidth, jittershape, jittersize, jittercol, names , log) {
  # Initiate the plot and add above LOD scatter points
  ymin <- if(lod <= min(pl$yi)) {
            lod;
          }
          else {
            lod-((1-max(summary(factor(pl$xi[pl$yi<lod]))) %% stack/stack) + 
                max(summary(factor(pl$xi[pl$yi<lod]))) / stack) *
                (max(pl$yi) - lod) / (1/ratio);
          };

  plot(pl$xi[pl$yi>=lod] + runif(length(pl$xi[pl$yi>=lod]), -jitterwidth, jitterwidth), pl$yi[pl$yi>=lod],
       xlim = c(0.25, vnum + 0.75), ylim = c(ymin, max(pl$yi)),
       xaxt = "n", yaxt = "n",
       main = main, xlab = xlab, ylab = ylab,
       pch = jittershape, cex = jittersize, col = jittercol
       );

  # Custom axis
  ticks <- seq(from = lod, to = max(pl$yi), length.out = axis);
  ticks[2:axis] <- round(ticks[2:axis], 2);
  axis(side = 2, at = ticks, if(log=="y") { labels = round(exp(ticks),2);});

  # Labels for each variable
  mtext(names, side = 1, at = seq(1:vnum));

  # LOD dashed line
  segments(0, lod, vnum + 1, lod, lty="dashed");
}
