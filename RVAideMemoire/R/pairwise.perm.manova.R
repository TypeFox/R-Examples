pairwise.perm.manova <- function(resp,fact,test=c("Pillai","Wilks","Hotelling-Lawley","Roy","Spherical"),
  nperm=999,p.method="fdr") {
  call <- match.call()
  dname <- paste0(deparse(call$resp)," by ",deparse(substitute(fact)),"\n",nperm," permutations")
  if ("dist" %in% class(resp)) {
    fun.p <- function(i,j) {
	fact2 <- droplevels(fact[as.numeric(fact)%in%c(i,j)])
	resp2 <- as.matrix(resp)
	rows <- which(fact%in%levels(fact2))
	resp2 <- as.dist(resp2[rows,rows])
	vegan::adonis(resp2~fact2,permutations=nperm)$aov.tab[1,"Pr(>F)"]
    }
    multcomp <- pairwise.table(fun.p,levels(fact),p.adjust.method=p.method)
    method <- "permutational MANOVAs on a distance matrix"
  } else {
    if (nrow(resp)!=length(fact)) {
	stop(paste("'",deparse(substitute(resp)),"' and '",deparse(substitute(fact)),
	  "' lengths differ",sep=""))
    }
    test <- match.arg(test)
    if (!is.matrix(resp)) {resp <- as.matrix(resp)}
    if (!is.factor(fact)) {fact <- factor(fact)}
    fun.p <- function(i,j) {
	resp2 <- resp[as.numeric(fact)%in%c(i,j),]
	fact2 <- droplevels(fact[as.numeric(fact)%in%c(i,j)])
	perm.manova(resp2,fact2,test=test,nperm=nperm)
    }
    multcomp <- pairwise.table(fun.p,levels(fact),p.adjust.method=p.method)
    method <- paste0("permutational MANOVAs (test: ",test,")")
  }
  result <- list(method=method,data.name=dname,p.value=multcomp,p.adjust.method=p.method)
  class(result) <- "pairwise.htest"
  return(result)
}

perm.manova <- function(resp,fact,test,nperm) {
  manova.ref <- anova(lm(resp~fact),test=test)
  stat <- ifelse(test!="Spherical","approx F","F")
  F.ref <- manova.ref[2,stat]
  F.perm <- numeric(nperm+1)
  F.perm[1] <- F.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    manova.perm <- anova(lm(resp[sample(1:nrow(resp)),]~fact),test=test)
    F.perm[i+1] <- manova.perm[2,stat]
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue <- length(which((F.perm+.Machine$double.eps/2) >= F.ref))/(nperm+1)
  return(pvalue)
}

