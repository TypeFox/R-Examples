#
# Limitplot functions -- plotBox
#

plotBox <- function(pl, lod, vnum, CI, blod) {
  # Add box plot(s)
  for(i in 1:vnum) {
    vect <- c(pl$yi[pl$yi>=lod & pl$xi==i], seq(from = lod * blod, to = lod * blod, length.out = length(pl$yi[pl$yi<lod & pl$xi==i])));
    avg <- mean(vect);
    stdev <- sd(vect) / sqrt(length(pl$yi[pl$xi==i]));
    rect(i-0.25, 
         max(lod, qnorm(((1-CI/100) / 2), avg, stdev)), 
         i+0.25, 
         max(lod, qnorm((1-(1-CI/100)/2), avg, stdev))
    );
  
    if(lod<=avg) { segments(i-0.25, avg, i+0.25, avg); }
  }
}
