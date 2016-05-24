##' function to make a simple set of graphs for a data set
##'
##' @param x data
##' @return NULL
##' 
##' @export
"simple.eda" <-
  function(x) {
    ### create a simple function to explore data
    op <- par(no.readonly = TRUE); # save old parameters
    par(mfrow=c(1,3))
    hist(x);rug(x)
    boxplot(x);rug(x,side=2);title("boxplot")
    qqnorm(x);qqline(x)
    par(op)
  }
