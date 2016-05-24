"plot.mrpp" <- 
function(x, y, ...) {
   with(x, {
   fig.dist <- hist(boot.deltas, xlim=c(range(delta, boot.deltas), ...), 
                  main="Are groups significantly different?", xlab=call)
   abline(v=delta); 
   text(delta, diff(range(fig.dist$counts))/2, adj = -0.5,
     expression(bold(delta)), cex=1.5 )  }
)
}