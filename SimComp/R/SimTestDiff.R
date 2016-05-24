SimTestDiff <-
function(data,grp,resp=NULL,type="Dunnett",base=1,ContrastMat=NULL,
                        alternative="two.sided",Margin=NULL,covar.equal=FALSE) {


alternative <- match.arg(alternative, choices=c("two.sided", "less", "greater"))
if (length(grp) > 1) {
  stop("Specify only one grouping variable")
}
tr.names <- levels(data[,grp])
ntr <- length(levels(data[,grp]))                                   # number of treatments
trlist <- split(x=data,f=data[,grp])                                # data list splitted into the different treatments
ssvec <- numeric(ntr)                                               # sample sizes for the treatments
for (i in 1:ntr) {
  ssvec[i] <- nrow(trlist[[i]])
}

if (is.null(resp)) {
  resp <- names(data[,-which(names(data)%in%grp),drop=FALSE])
} 
nep <- length(resp)                                                 # number of endpoints

if (covar.equal==TRUE) {
  if ( sum(ssvec-1)<nep ) {
    stop("Not enough observations or too many endpoints to be analyzed")
  }
} else {
  if ( any((ssvec-1)<nep) ) {
    stop("Not enough observations or too many endpoints to be analyzed")
  }
}

for (i in 1:ntr) {
  trlist[[i]] <- as.matrix(trlist[[i]][,resp])                      # data list without factor
}

for (i in 1:ntr) {
  if (is.numeric(trlist[[i]])==FALSE) {
    stop("Response variables must be numeric")
  }
  if (any(is.na(trlist[[i]]))) {
    stop("Unequal sample sizes for the endpoints; missing values")
  }
}

if (!is.null(ContrastMat)) {
  if (ncol(ContrastMat)!=ntr) {
    stop("Number of columns of ContrastMat and number of groups must be equal")
  }
  C0 <- apply(X=ContrastMat, MARGIN=1, function(x) {
    all(x==0)
  })
  if (any(C0)) {
    cat("Warning: At least one row of ContrastMat is a vector with all components", "\n",
        "equal to zero", "\n")
  }
  Cmat <- ContrastMat
  type <- "User defined"
  if (is.null(rownames(Cmat))) {
    rownames(Cmat) <- paste("C", 1:nrow(Cmat), sep="")
  }
} else {
  type <- match.arg(type, choices=c("Dunnett", "Tukey", "Sequen", "AVE", "GrandMean", "Changepoint", 
    "Marcus", "McDermott", "Williams", "UmbrellaWilliams"))
  names(ssvec) <- tr.names
  Cmat <- contrMat(n=ssvec, type=type, base=base)
}
comp.names <- rownames(Cmat)

if (is.null(Margin)) {
  Margin <- matrix(rep(0,nrow(Cmat)*nep),ncol=nep)
}
if (is.numeric(Margin)) {
  if (length(Margin)==1) {
    Margin <- matrix(rep(Margin,nrow(Cmat)*nep),ncol=nep)
  }
  if (is.vector(Margin)) {
    if (length(Margin)!=nep) {
      stop("Margin must be a single numeric value, or a numeric vector with length equal to the number of endpoints,", "\n",
           "or a matrix with rows equal to the number of contrasts and columns equal to the number of endpoints")
    } else {
      Margin <- matrix(rep(Margin,nrow(Cmat)),ncol=nep,byrow=TRUE)
    }
  }
  if (is.matrix(Margin)) {
    if (nrow(Margin)!=nrow(Cmat)) {
      stop("Number of rows of Margin and number of contrasts must be equal")
    }
    if (ncol(Margin)!=nep) {
      stop("Number of columns of Margin and number of endpoints must be equal")
    }
  }
}

if (covar.equal==TRUE) {
  out <- SimTestDiffHom(trlist=trlist, grp=grp, ntr=ntr, nep=nep, ssvec=ssvec, Cmat=Cmat,
                        alternative=alternative, Margin=Margin)
} else {
  out <- SimTestDiffHet(trlist=trlist, grp=grp, ntr=ntr, nep=nep, ssvec=ssvec, Cmat=Cmat,
                        alternative=alternative, Margin=Margin)
}
out$type <- type
out$test.class <- "differences"
out$covar.equal <- covar.equal
out$comp.names <- comp.names
out$resp <- resp
colnames(out$estimate) <- resp
if (type=="User defined") {
  rownames(out$estimate) <- comp.names
  colnames(out$Cmat) <- tr.names
}
rownames(out$statistic) <- comp.names; colnames(out$statistic) <- resp
rownames(out$p.val.raw) <- comp.names; colnames(out$p.val.raw) <- resp
rownames(out$p.val.adj) <- comp.names; colnames(out$p.val.adj) <- resp
rownames(out$Margin) <- comp.names; colnames(out$Margin) <- resp
if (covar.equal==FALSE) {
  names(out$degr.fr) <- comp.names
  names(out$CovMatDat) <- names(out$CorrMatDat) <- tr.names
}
class(out) <- "SimTest"
return(out)


}
