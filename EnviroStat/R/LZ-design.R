
# This R function calculates the log |determinant| off all sub-covariance
# matrices of size (k x k) from a covariance matrix.

ldet.eval = function(covmat,k,all = FALSE)
# Input:
#   covmat: a covariance matrix (ie. non-negative definite, 
#           square and symmetric)
#   k: dimension of sub-covariance matrices considered
#  Optional: if True, returns all combinations with corresponding log|det|
#    Note - This option may need additionally a large amount of memory and
#            so may not work for a large number of combinations!!
# Output: 
#   coord.sel: The k coordinates having the largest log|det|
#   log.det  : The log|det| of the submatrix corresponding the coord.sel
#   all.comb : Null if all = False
#              all combinations and their log|det| if all = True
 { 
   n= dim(covmat)[1]
   if (k > n) stop("k is larger than the dimension of the matrix")
   if (sum(covmat == t(covmat)) < n^2) {
    warning("This R function only works for non-negative definite and symmetric matrix")
    warning("    The data matrix provided may not be symmetric")
         }
   tmp = eigen(covmat)$values
   if ((mode(tmp) == "numeric")&(min(tmp) < 0)) 
     warning("  The data matrix may not be non-negative definite")
   covmatv = c(covmat)
 if (!all) {
   all.comb = NULL
   icomb = rep(0, k)
   selcomb = rep(0,k)   
   numcomb =0
   ldet = 0
   p = rep(0,n +2)
   b = rep(0,n) 
   tmp = .C("search_ent",
   as.double(covmatv),
   as.integer(n),
   as.integer(k),
   as.integer(icomb),
   numcomb = as.integer(numcomb),
   selcomb = as.integer(selcomb),
   as.integer(p),
   as.integer(b),
   ldet = as.double(ldet))
   out =list(coord.sel=tmp$selcomb,log.det=tmp$ldet,all.comb = all.comb)
  }
 else
 {
   ncomb = choose(n,k)
   comb = rep(0, k*ncomb)
   p = rep(0,n +2)
   b = rep(0,n) 
   comball = .C("combgen",
                as.integer(n),
                as.integer(k),
                comb = as.integer(comb),
                as.integer(p),
                as.integer(b))$comb

   ldet = rep(0,ncomb)
   log.det = .C("eval_ent",
               as.double(covmatv),
               as.integer(n),
               as.integer(k),
               as.integer(comball),
               as.integer(ncomb),
               ldet = as.double(ldet))$ldet
   tmp = cbind(matrix(comball,byrow=T,ncol=k),log.det)
   ord = order(log.det[1:ncomb],decreasing =T)
   tmp = tmp[ord,]
   out=list(coord.sel=tmp[1,-(k+1)],log.det=tmp[1,k+1],all.comb = tmp)
  }
  return(out)
 }


# This R function lists all combinations of size 'k' from 
# the set {1,.., n}.
comb.all = function(n,k)
 {
   ncomb = choose(n,k)
   comb = rep(0, k*ncomb)
   p = rep(0,n +2)
   b = rep(0,n) 
   tmp = .C("combgen",
   as.integer(n),
   as.integer(k),
   comb = as.integer(comb),
   as.integer(p),
   as.integer(b))$comb
   all.comb = matrix(tmp,byrow=T,ncol=k)
}

