`histpp` <-
function(x, xname="", utitle="") {
  hist(x,main=utitle,xlab="", col="lightgrey", ylab=xname)
  rug(x, side=1, ticksize=0.03,  col="red")
  rug(x, side=1, ticksize=0.01,  col="white")
  rug(x, side=1, ticksize=0,  col="black")
}

