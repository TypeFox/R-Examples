fs.leave.k.out.sample <-
function(r.idx,dt,classlabel,
           k=2,n,
           method.agreement=c("euclidean","manhattan","Rand","cRand","NMI",
             "KP","angle","diag","Jaccard","FM","gdbr"),
           method.dist=c("pearson","kendall","spearman","standardizedEuclid",
             "euclidean","euclidean","pearson.u","kendall.u","spearman.u"),
           leave.by=c("count.class","flat","percent.class"),leave.k=1,
           return.value=c("mean","sd"))
  ## leave k out algorithm with exactly every possible combination for classlabel
  {

    method.agreement <- match.arg(method.agreement)
    method.dist <- match.arg(method.dist)
    return.value <- match.arg(return.value)
    if(is.null(r.idx)) r.idx <- seq(nrow(dt))
    leave.by <- match.arg(leave.by)
    if(identical(leave.by,"flat")) {
      nu.idx <- matrix(sample(seq(classlabel),
                              (length(classlabel)-leave.k)*n,
                              replace=T),
                       nrow=n)
    }
    else if(identical(leave.by,"count.class")) {
      nu.idx <- sample.classlabel.idx(classlabel,leave.k=leave.k,leave.k.mode="count",n=n)
    }
    else if(identical(leave.by,"percent.class")) {
      nu.idx <- sample.classlabel.idx(classlabel,
                                      leave.k=leave.k,
                                      leave.k.mode="percent",
                                      n=n)
    }
    else stop("Unsupported leave.by value!")
    agreement.measure <- apply(nu.idx,
                               1,
                               fs.agreement.part,
                               r.idx=r.idx,
                               dt=dt,
                               classlabel=classlabel,
                               k=k,
                               method.agreement=method.agreement,
                               method.dist=method.dist)
    if(identical(return.value,"mean")) return(mean(agreement.measure))
    else if(identical(return.value,"sd")) return(sd(agreement.measure))
    else stop("Unsupported return values")
  }
