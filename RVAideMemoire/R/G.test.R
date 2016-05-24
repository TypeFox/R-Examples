G.test <- 
function(x,p=rep(1/length(x),length(x))) {
  if (any(x==0)) {warning("G test should not be used with 0 values")}
  call <- match.call()
  data.name <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
  if (is.table(x) | is.matrix(x)) {
    chi <- suppressWarnings(chisq.test(x))
    method <- "G-test"
  } else if (is.vector(x)) {
    chi <- suppressWarnings(chisq.test(x,p=p))
    method <- "G-test for given probabilities"
  }
  theo <- chi$expected
  x2 <- x[x>0]
  theo2 <- theo[x>0]
  G <- 2*sum(x2*log(x2/theo2))
  names(G) <- "G"
  ddl <- chi$parameter
  p <- pchisq(G,ddl,lower.tail=FALSE)
  result <- list(method=method,statistic=G,parameter=ddl,p.value=p,
    data.name=data.name,observed=x,expected=theo)
  class(result) <- "htest"
  return(result)
}
