print.MVA.test.mult <- function(x,digits=4,...) {
  cat("\n")
  cat(strwrap(x$method,prefix="\t"),sep="\n")
  cat("\n")
  cat("data: ",x$data.name,"\n\n")
  tab <- data.frame(Response=names(x$p.value),Q2=x$statistic,P=x$p.value)
  print(tab,digits=digits,row.names=FALSE)
  cat(paste("\nP value adjustment method: ",x$p.adjust.method,"\n",sep=""))
  invisible(x)
}

MVA.test <- function(X,Y,cmv=FALSE,ncomp=5,kout=7,kinn=8,model=c("PLSR","CPPLS","PLS-DA","PPLS-DA","LDA",
  "QDA","PLS-DA/LDA","PLS-DA/QDA","PPLS-DA/LDA","PPLS-DA/QDA"),Q2diff=0.05,lower=0.5,upper=0.5,Y.add=NULL,
  weights=rep(1,nrow(X)),set.prior=FALSE,crit.DA=c("plug-in","predictive","debiased"),p.method="fdr",nperm=999,
  ...) {
  model <- match.arg(model)
  dname <- paste0(deparse(substitute(X))," and ",deparse(substitute(Y)),"\nModel: ",model,
    "\n",ncomp," components",ifelse(cmv," maximum",""),"\n",nperm," permutations")
  if (!cmv) {
    testname <- "Permutational test based on cross-validation"
    cv.ref <- MVA.cv(X,Y,repet=20,k=kout,ncomp=ncomp,model=model,lower=lower,upper=upper,Y.add=Y.add,
	weights=weights,set.prior=set.prior,crit.DA=crit.DA,...)
  } else {
    testname <- "Permutational test based on cross model validation"
    cv.ref <- MVA.cmv(X,Y,repet=20,kout=kout,kinn=kinn,ncomp=ncomp,model=model,crit.inn=ifelse(is.factor(Y),"NMC","Q2"),
	Q2diff=Q2diff,lower=lower,upper=upper,Y.add=Y.add,weights=weights,set.prior=set.prior,crit.DA=crit.DA,...)
  }
  if(cv.ref$type=="quant") {
    ref <- colMeans(cv.ref$Q2)
    names(ref) <- if (ncol(as.data.frame(Y))==1) {"Q2"} else {paste("Q2",colnames(Y),sep=".")}
  } else {
    ref <- mean(as.vector(cv.ref$NMC))
    names(ref) <- "NMC"
  }
  if (cv.ref$type=="quant" & ncol(as.data.frame(Y))>1) {
    stat.perm <- matrix(0,nrow=nperm+1,ncol=ncol(Y))
    stat.perm[1,] <- ref
  } else {
    stat.perm <- numeric(nperm+1)
    stat.perm[1] <- ref
  }
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for(i in 1:nperm) {
    setTxtProgressBar(pb,round(i*100/nperm,0))
    if (!cmv) {
	if (cv.ref$type=="quant" & ncol(as.data.frame(Y))>1) {
	  cv.perm <- MVA.cv(X,Y[sample(1:nrow(Y)),],repet=1,k=kout,ncomp=ncomp,model=model,lower=lower,upper=upper,
	    Y.add=Y.add,weights=weights,set.prior=set.prior,crit.DA=crit.DA,...)
	} else {
	  cv.perm <- suppressWarnings(MVA.cv(X,sample(Y),repet=1,k=kout,ncomp=ncomp,model=model,lower=lower,
	    upper=upper,Y.add=Y.add,weights=weights,set.prior=set.prior,crit.DA=crit.DA,...))
	}
    } else {
	if (cv.ref$type=="quant" & ncol(as.data.frame(Y))>1) {
	  cv.perm <- MVA.cmv(X,Y[sample(1:nrow(Y)),],repet=1,kout=kout,kinn=kinn,ncomp=ncomp,model=model,
	    crit.inn=ifelse(is.factor(Y),"NMC","Q2"),Q2diff=Q2diff,lower=lower,upper=upper,
	    Y.add=Y.add,weights=weights,set.prior=set.prior,crit.DA=crit.DA,...)
	} else {
	  cv.perm <- suppressWarnings(MVA.cmv(X,sample(Y),repet=1,kout=kout,kinn=kinn,ncomp=ncomp,model=model,
	    crit.inn=ifelse(is.factor(Y),"NMC","Q2"),Q2diff=Q2diff,lower=lower,upper=upper,Y.add=Y.add,
	    weights=weights,set.prior=set.prior,crit.DA=crit.DA,...))
	}
    }
    if (cv.ref$type=="quant" & ncol(as.data.frame(Y))>1) {
	stat.perm[i+1,] <- cv.perm$Q2
    } else {
	stat.perm[i+1] <- if(cv.ref$type=="quant") {as.vector(cv.perm$Q2)} else {as.vector(cv.perm$NMC)}
    }
  }
  if (cv.ref$type=="quant") {
    if (ncol(as.data.frame(Y))==1) {
	pvalue <- length(which((stat.perm+.Machine$double.eps/2) >= ref))/(nperm+1)
    } else {
	pvalue <- numeric(ncol(Y))
	for (i in 1:ncol(Y)) {pvalue[i] <- length(which((stat.perm[,i]+.Machine$double.eps/2) >= ref[i]))/(nperm+1)}
	pvalue <- p.adjust(pvalue,method=p.method)
	names(pvalue) <- colnames(Y)
    }
  } else {
    pvalue <- length(which((stat.perm-.Machine$double.eps/2) <= ref))/(nperm+1)
  }
  result <- list(method=testname,data.name=dname,statistic=ref,permutations=nperm,p.value=pvalue,
    p.adjust.method=p.method)
  class(result) <- ifelse(cv.ref$type=="quant" & ncol(as.data.frame(Y))>1,"MVA.test.mult","htest")
  return(result)
}


