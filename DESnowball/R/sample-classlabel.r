sample.classlabel.idx <-
function(classlabel,leave.k,n,leave.k.mode)
  ##
{
    alist <- sample.classlabel.idx.aux2(classlabel=classlabel,
					leave.k=leave.k,
					n=n,
					leave.k.mode=leave.k.mode)
    for (i in seq(alist)) {
	if(i == 1) ret <- alist[[i]]
	else ret <- cbind(ret,alist[[i]])
    }
    ret
}

sample.classlabel.idx.aux1 <-
function(x,leave.size,n,...)
  ## sample into a matrix format
{
    matrix(sample(x,(length(x)-leave.size)*n,...),nrow=n)
}

sample.classlabel.idx.aux2 <-
function(classlabel,leave.k,n,leave.k.mode=c('count','percent'))
  ## n is the sample size
{
    classlabel <- factor(classlabel)
    leave.k.mode <- match.arg(leave.k.mode)
    if(length(leave.k)!=nlevels(classlabel) & length(leave.k) != 1)
	stop("length of leave.k needs to be the same as nlevels of classlabel or 1")
    else if (length(leave.k)==1 & identical(leave.k.mode,"count")) {
	ret <- tapply(seq(classlabel),
		      classlabel,
		      sample.classlabel.idx.aux1,
		      leave.size=leave.k,
		      n=n,
		      replace=T)
    }
    else {
	ret <- list()
	for (i in seq(nlevels(classlabel))) {
	    if(identical(leave.k.mode,"percent")) {
		if(length(leave.k) == 1)
		    m <- floor(sum(classlabel==levels(classlabel)[i])*leave.k)
		else m <- floor(sum(classlabel==levels(classlabel)[i])*leave.k[i])
	    }
	    else m <- leave.k[i]
	    ret[[i]] <- sample.classlabel.idx.aux1(which(classlabel==levels(classlabel)[i]),
						   leave.size=m,n=n,replace=T)
	}
    }
    names(ret) <- levels(classlabel)
    ret
}
