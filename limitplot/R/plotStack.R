#
# Limitplot functions -- plotStack
#

plotStack <- function(pl, lod, vnum, stack, ratio, shape, size, col) {
  # Stack the point(s) below the LOD
  if(lod>min(pl$yi)) {
    for(i in 1:vnum) {
      xp<-rep(seq(-0.2, 0.2, length.out = stack), len = length(pl$xi[pl$xi==i & pl$yi<lod])) + i;
      yp<-rep(seq(1:(1-max(summary(factor(pl$xi[pl$yi<lod]))) %% stack / stack) + 
              max(summary(factor(pl$xi[pl$yi<lod]))) / stack), each = stack, len = length(pl$xi[pl$yi<lod & pl$xi==i]));
      points(xp, lod - yp * (max(pl$yi) - lod) / (1/ratio), pch = shape, cex = size, col = col);
    }
  }
}


