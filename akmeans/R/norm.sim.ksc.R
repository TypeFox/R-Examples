norm.sim.ksc <-
function(A,k,init.cen=NULL,init.mem=NULL,iter.max=100){
## assume A is already normalized by square root of sum=> every row has L2norm as 1.
  n = nrow(A); p = ncol(A)
## centeroids init
  if (!is.null(init.cen)) {
    cur.cen = init.cen
  } else if (!is.null(init.mem)) {
    cur.cen = norm.sim.ksc.center.update(init.mem,A,k)
  } else {cur.cen = A[sample(n,k),]}

## assignment init
  cur.mem = apply(A%*%t(cur.cen),1,which.max)
  for (i in 1:iter.max){
    prev.mem = cur.mem
    cur.cen = norm.sim.ksc.center.update(cur.mem, A,k)
    cur.mem = apply(A%*%t(cur.cen),1,which.max)
    if (sum(abs(prev.mem-cur.mem))==0) break ## convergence check
  }
  cur.cen = norm.sim.ksc.center.update(cur.mem, A,k)
  size = apply(matrix(1:k,ncol=1),1,function(i){sum(cur.mem==i)})
  list(cluster = cur.mem, centers = cur.cen, size = size)
}

