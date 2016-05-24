SimCiDiff <-
function(data,grp,resp=NULL,type="Dunnett",base=1,ContrastMat=NULL,
                       alternative="two.sided",covar.equal=FALSE,conf.level=0.95) {


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

if (length(conf.level)>1 || !is.numeric(conf.level) || conf.level>=1 || conf.level<=0 ) {
  stop("conf.level must be a single numeric value between 0 and 1")
}

if (covar.equal==TRUE) {
  out <- SimCiDiffHom(trlist=trlist, grp=grp, ntr=ntr, nep=nep, ssvec=ssvec, Cmat=Cmat,
                       alternative=alternative, conf.level=conf.level)
} else {
  out <- SimCiDiffHet(trlist=trlist, grp=grp, ntr=ntr, nep=nep, ssvec=ssvec, Cmat=Cmat,
                       alternative=alternative, conf.level=conf.level)
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
rownames(out$lower) <- comp.names; colnames(out$lower) <- resp
rownames(out$upper) <- comp.names; colnames(out$upper) <- resp
if (covar.equal==FALSE) {
  names(out$degr.fr) <- comp.names
  names(out$CovMatDat) <- names(out$CorrMatDat) <- tr.names
}
class(out) <- "SimCi"
return(out)


}
