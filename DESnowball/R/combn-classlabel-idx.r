combn.classlabel.idx <-
function(classlabel,leave.k,leave.k.mode)
  ## combn on classlabel in different ways
{
    nu.idx.byclass <- combn.classlabel.idx.aux1(classlabel=classlabel,
						leave.k=leave.k,
						leave.k.mode=leave.k.mode)
    expand.idx <- expand.grid(lapply(nu.idx.byclass,seq.ncol))
    extract.classlabel <- function(eachRow.in.expand.idx,
				   nu.idx.byclass,classlabel)
	seq(classlabel)[-unlist(combn.classlabel.idx.aux2(nu.idx.byclass,
							  eachRow.in.expand.idx))]    
    ret <- t(apply(expand.idx,1,extract.classlabel,nu.idx.byclass,classlabel))
    ret
}

combn.classlabel.idx.aux1 <-
function(classlabel,leave.k,leave.k.mode=c("count","percent"))
  ## tapply with combn but with different count or percentage for each level of classlabel
{
    classlabel <- factor(classlabel)
    leave.k.mode <- match.arg(leave.k.mode)
    if(length(leave.k)!=nlevels(classlabel) & length(leave.k) != 1)
	stop("length of leave.k needs to be the same as nlevels of classlabel or 1")
    else if (length(leave.k)==1 & identical(leave.k.mode,"count")) {
	ret <- tapply(seq(classlabel),classlabel,combn,m=leave.k)
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
	    ret[[i]] <- combn(which(classlabel==levels(classlabel)[i]),m)
	}
    }
    names(ret) <- levels(classlabel)
    ret
}

combn.classlabel.idx.aux2 <-
function(x,e,...)
  ## lapply to the list x, each elemen in x is applied with one element in e,
  ##  so e and x need to have the same length
{
    ret <- list()
    for (i in seq(x)) {
	ret[[i]] <- x[[i]][,e[i],...]
    }
    ret
}

