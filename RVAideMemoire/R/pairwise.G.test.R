pairwise.G.test <-
function(x,p.method="fdr") {
  fun.p <- function(i,j) {
    y <- x[c(i,j),]
    G.test(y)$p.value
  }
  level.names <- if (!is.null(rownames(x))) {
    rownames(x)
  } else {
    LETTERS[1:nrow(x)]
  }
  tab.p <- pairwise.table(fun.p,level.names,p.adjust.method=p.method)
  call <- match.call()
  dname.x <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
  result <- list(method="G-tests",data.name=dname.x,p.adjust.method=p.method,p.value=tab.p)
  class(result) <- "pairwise.htest"
  return(result)
}
