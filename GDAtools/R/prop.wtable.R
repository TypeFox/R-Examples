prop.wtable <- function(var1,var2=NULL,w=rep.int(1,length(var1)),dir=0,digits=1,mar=TRUE,na=TRUE) {
  t <- wtable(var1,var2,w=w,digits=10,mar=TRUE,na=na)
  if(is.null(var2)) {
    wtab <- 100*2*t/sum(t)
    rownames(wtab) <- rownames(t)
    if(mar==FALSE) wtab <- as.matrix(wtab[-length(wtab),])
  } else {
    if(dir==0) wtab <- 100*4*t/sum(t)
    if(dir==1) wtab <- apply(t,2,function(x) 100*2*x/rowSums(t))
    if(dir==2) wtab <- t(apply(t,1,function(x) 100*2*x/colSums(t)))
    dimnames(wtab) <- dimnames(t)
    if(mar==FALSE) wtab <- wtab[-nrow(wtab),-ncol(wtab)]
    }
  wtab <- round(wtab,digits)
  return(wtab)
  }
