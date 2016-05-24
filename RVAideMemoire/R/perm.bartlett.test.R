perm.bartlett.test <-
function(formula,data,nperm=999) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  m <- match.call()
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$nperm <- NULL
  mf <- eval(m,parent.frame())
  dname <- paste(paste(names(mf)[1],paste(names(mf)[2:ncol(mf)],collapse=":"),sep=" by "),"\n",nperm," permutations",sep="")
  resp <- mf[,1]
  fact <- interaction(mf[,2:ncol(mf)],sep=":")
  K.ref <- bartlett.test(resp~fact)$statistic
  K.perm <- numeric(nperm+1)
  K.perm[1] <- K.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for(i in 1:nperm) {
    K.perm[i+1] <- bartlett.test(sample(resp)~fact)$statistic
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue <- length(which((K.perm+.Machine$double.eps/2) >= K.ref))/(nperm+1)
  result <- list(method="Permutational Bartlett test of homogeneity of variances",data.name=dname,statistic=K.ref,
    permutations=nperm,p.value=pvalue)
  class(result) <- "htest"
  return(result)
}
