simpute.svd <-
  function (x, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=FALSE,warm.start=NULL,...) 
{
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  xnas <- is.na(x)

  nz=m*n-sum(xnas)
  xfill <- x
  xfill[xnas] <- 0
  if(!is.null(warm.start)){
###must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    D=warm.start$d
    nzD=sum(D>0)
    JD=min(nzD,J)
    U=warm.start$u[,seq(JD),drop=FALSE]
    V=warm.start$v[,seq(JD),drop=FALSE]
    D=D[seq(JD)]
    xhat=U%*%(D*t(V))
    xfill[xnas] <- xhat[xnas]
    }
  svd.xfill=svd(xfill)
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    svd.old=svd.xfill
    d=svd.xfill$d
    d=pmax(d-lambda,0)
    xhat <- svd.xfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.xfill$v[,seq(J)]))
    xfill[xnas] <- xhat[xnas]
    svd.xfill=svd(xfill)
    ratio=Frob(svd.old$u[, seq(J)],d[seq(J)],svd.old$v[, seq(J)],
      svd.xfill$u[, seq(J)],pmax(svd.xfill$d-lambda,0)[seq(J)],svd.xfill$v[, seq(J)])
    if(trace.it){
      obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(d))/nz
      cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
    }
  }
  d=pmax(svd.xfill$d[seq(J)]-lambda,0)
  J=min(sum(d>0)+1,J)
  svd.xfill=list(u=svd.xfill$u[, seq(J)], d=d[seq(J)], v=svd.xfill$v[,seq(J)])
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  attributes(svd.xfill)=c(attributes(svd.xfill),list(lambda=lambda,call=this.call),biats)
  svd.xfill
}
