plotcoverage<-function(x, y=NULL, by=0.05, type="stat", max.order=NULL, group = NULL,  alpha = 0.05, add=FALSE, xlab=expression(A[t]),...) {
  match.arg(type,c("lowerCI","upperCI","stat"))
  Atseq <- seq(0,1.0,by=by)
  covA<-numeric(length(Atseq))
  if(inherits(x,"indicators")) {     
    num.order <- rowSums(x$C)
    sel2 <- rep(TRUE,length(num.order))
    if(!is.null(max.order)) sel2[num.order>max.order] <- FALSE
    for(i in 1:length(Atseq)) {
      covA[i] = coverage(x, selection = sel2, At = Atseq[i], type=type)
    }
  }
  else if(inherits(x,"data.frame")) {
    if(is.null(y)) stop("You must supply an object of class 'multipatt' for parameter 'y'")
    if(!inherits(y,"multipatt")) stop("Wrong class for 'y'. Should be an object of class 'multipatt'")
    if(is.null(group)) stop("Please, provide a value for parameter 'group': a name or number for a site group or site group combination")     
    if(is.numeric(group)) combnumber = group
    else combnumber = which(colnames(y$comb)==group)
    if(length(combnumber)==0) stop("Wrong group name")
    for(i in 1:length(Atseq)) {
      covA[i] = coverage(x, y, At = Atseq[i], alpha = alpha)[combnumber]
    }
  }
  if(!add) {
    plot(Atseq,covA*100, type="s", ylim=c(0,100), axes=FALSE, xlab=xlab, ylab="Coverage (%)", ...)  
    axis(1, pos=0)
    axis(2, pos=0)     
  } else lines(Atseq,covA*100, type="s", ...)
}
