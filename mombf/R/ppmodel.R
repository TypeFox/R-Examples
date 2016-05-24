###
### ppmodel.R
###

ppmodel <- function(nlpfit) {
  postmodel <- apply(nlpfit$postModel==1,1,function(z) paste(which(z),collapse=','))
  ans <- table(postmodel)
  ans <- ans[order(ans,decreasing=TRUE)]
  nvars <- sapply(strsplit(names(ans),split=','),length)
  ans <- data.frame(selectVars=names(ans),nvars=nvars,count=ans,posprob=ans/sum(ans))
  rownames(ans) <- NULL
  return(ans)
}

