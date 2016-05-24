"plot.dsl" <- 
function(x, y, ...) {
   with(x, {
   fig.dist <- hist(perms, xlim=c(range(slope.diff, perms), ...), 
                  main="Is difference in slope significant?", xlab=call)
   abline(v=slope.diff); 
   text(slope.diff, diff(range(fig.dist$counts))/2, adj = -0.5,
     expression(bold(ds0)), cex=1.5 )  }
)
}