pairwise.mood.medtest <- function(resp,fact,exact=NULL,p.method="fdr") {
  if (length(resp)!=length(fact)) {
    stop(paste("'",deparse(substitute(resp)),"' and '",deparse(substitute(fact)),
	"' lengths differ",sep=""))
  }
  if (!is.numeric(resp)) {resp <- as.numeric(as.character(resp))}
  if (!is.factor(fact)) {fact <- factor(fact)}
  data.name <- paste(deparse(substitute(resp)),"and",deparse(substitute(fact)))
  method <- "Mood's median tests"
  if(is.null(exact)) {exact <- length(resp)<200}
  fun.p <- function(i,j) {
    resp2 <- resp[as.numeric(fact)%in%c(i,j)]
    fact2 <- droplevels(fact[as.numeric(fact)%in%c(i,j)])
    mood.medtest(resp2~fact2,exact=exact)$p.value
  }
  multcomp <- pairwise.table(fun.p,levels(fact),p.adjust.method=p.method)
  result <- list(method=method,data.name=data.name,p.value=multcomp,p.adjust.method=p.method)
  class(result) <- "pairwise.htest"
  return(result)
}
