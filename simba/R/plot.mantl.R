"plot.mantl" <- 
function(x, y, ...) {
   with(x, {
   fig.dist <- hist(perms, xlim=c(range(statistic, perms), ...), 
                  main="Is correlation significant?", xlab=call)
   abline(v=statistic); 
   text(statistic, diff(range(fig.dist$counts))/2, adj = -0.5,
     expression(bold(ds0)), cex=1.5 )  }
)
}