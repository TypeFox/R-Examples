chisq.multcomp <-
function(x,p.method="fdr") {
  x <- sort(x)
  fun.p <- function(i,j) {
    xi <- x[i]
    xj <- x[j]
    suppressWarnings(chisq.test(c(xi, xj)))$p.value
  }
  tab.p <- pairwise.table(fun.p,as.character(x),p.adjust.method=p.method)
  call <- match.call()
  dname.x <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
  result <- list(method="chi-squared tests",data.name=dname.x,p.adjust.method=p.method,p.value=tab.p)
  class(result) <- "pairwise.htest"
  return(result)
}
