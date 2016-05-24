perm.cor.test <-
function(x,y,alternative=c("two.sided","less","greater"),nperm=999) {
  if (length(x)!=length(y)) {stop(paste("'",deparse(substitute(x)),"' and '",deparse(substitute(y)),"' lengths differ",sep=""))}
  if (!is.numeric(x)) {x <- as.numeric(as.character(x))}
  if (!is.numeric(y)) {y <- as.numeric(as.character(y))}
  if (length(alternative)>1) {alternative <- "two.sided"}
  data.name <- paste(deparse(substitute(x))," and ",deparse(substitute(y)),"\n",nperm," permutations",sep="")
  coeff <- cor.test(x,y,alternative=alternative)$estimate
  t.ref <- cor.test(x,y,alternative=alternative)$statistic
  t.perm <- numeric(nperm+1)
  t.perm[1] <- t.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for(i in 1:nperm) {
    setTxtProgressBar(pb,round(i*100/nperm,0))
    t.perm[i+1] <- cor.test(x,sample(y),alternative=alternative)$statistic
  }
  cat("\n")
  pvalue <- NULL
  if (alternative=="two.sided") {
    pvalue <- 2*min(length(which((t.perm-.Machine$double.eps/2) <= t.ref))/(nperm+1),length(which((t.perm+.Machine$double.eps/2) >= t.ref))/(nperm+1))
  }
  if (alternative=="less") {
    pvalue <- length(which((t.perm-.Machine$double.eps/2) <= t.ref))/(nperm+1)
    }
  if (alternative=="greater") {
    pvalue <- length(which((t.perm+.Machine$double.eps/2) >= t.ref))/(nperm+1)
  }
  null.value <- 0
  names(null.value) <- "correlation"
  result <- list(statistic=t.ref,permutations=nperm,p.value=pvalue,estimate=coeff,alternative=alternative,data.name=data.name,
    null.value=null.value,method="Pearson's product-moment correlation - Permutation test")
  class(result) <- "htest"
  return(result)
}
