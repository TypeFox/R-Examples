biScale<- function(x,maxit=20,thresh=1e-9,row.center=TRUE,row.scale=TRUE,col.center=TRUE,col.scale=TRUE,trace=FALSE){
  ### Function for doing both row and column centering and scaling
  ### Both generic, and methods for "sparseMatrix" and "Incomplete"
  mn=dim(x)
  m=mn[1];n=mn[2]
### Check row centering
     if(is.numeric(row.center)){
       if(length(row.center)==m){
         alpha=row.center
         row.center=FALSE
       }
       else stop("length of 'row.center' must equal the number of rows of 'x'")
     }
     else alpha=rep(0,m)
### Check column centering
     if(is.numeric(col.center)){
       if(length(col.center)==n){
         beta=col.center
         col.center=FALSE
       }
       else stop("length of 'col.center' must equal the number of columns of 'x'")
     }
     else beta=rep(0,n)
### Check row scaling
     if(is.numeric(row.scale)){
       if(length(row.scale)==m){
         if(any(row.scale<=0))stop("elements of 'row.scale' must be strictly positive")
         tau=row.scale
         row.scale=FALSE
       }
       else stop("length of 'row.scale' must equal the number of rows of 'x'")
     }
     else tau=rep(1,m)
### Check column scaling
     if(is.numeric(col.scale)){
       if(length(col.scale)==n){
         if(any(col.scale<=0))stop("elements of 'col.scale' must be strictly positive")
         gamma=col.scale
         col.scale=FALSE
       }
       else stop("length of 'col.scale' must equal the number of cols of 'x'")
     }
     else gamma=rep(1,n)

  out=centerScale(x,maxit,thresh,row.center,row.scale,col.center,col.scale,trace,m,n,alpha,beta,tau,gamma)
  critmat=attr(out,"critmat")
  if((nrow(critmat)==maxit)&&(critmat[maxit,2]>thresh))warning(paste("biScale not converged after",maxit,"iterations; try larger value of maxit and use trace=TRUE to monitor convergence"))
  out
}
