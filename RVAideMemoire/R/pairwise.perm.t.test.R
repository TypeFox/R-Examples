pairwise.perm.t.test <- function(resp,fact,p.method="fdr",paired=FALSE,
  alternative=c("two.sided","less","greater"),nperm=999) {
  if (length(resp)!=length(fact)) {
    stop(paste("'",deparse(substitute(resp)),"' and '",deparse(substitute(fact)),
	"' lengths differ",sep=""))
  }
  if (!is.numeric(resp)) {resp <- as.numeric(as.character(resp))}
  if (!is.factor(fact)) {fact <- factor(fact)}
  if (paired==TRUE & any(diff(tapply(resp,fact,length)))!=0) {
    stop(paste("all levels of ",deparse(substitute(fact))," must have the same length",sep=""))
  }
  if (length(alternative)>1) {alternative <- "two.sided"}
  data.name <- paste(deparse(substitute(resp))," and ",deparse(substitute(fact)),"\n",nperm," permutations",sep="")
  method <- if (!paired) {"permutational t tests"} else {"permutational paired t tests"}
  fun.p <- function(i,j) {
    resp2 <- resp[as.numeric(fact)%in%c(i,j)]
    fact2 <- droplevels(fact[as.numeric(fact)%in%c(i,j)])
    perm.t.test(resp2~fact2,alternative=alternative,paired=paired,nperm=nperm)$p.value
  }
  multcomp <- pairwise.table(fun.p,levels(fact),p.adjust.method=p.method)
  result <- list(method=method,data.name=data.name,p.value=multcomp,p.adjust.method=p.method,
    permutations=nperm)
  class(result) <- "pairwise.htest"
  return(result)
}
