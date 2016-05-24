is.matrix_ergodic <-
function(A,digits=12,return.eigvec=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
leigs<-eigen(t(A))
lmax<-which.max(Re(leigs$values))
v<-leigs$vectors[,lmax]
Rev<-abs(Re(v))
Rev<-round(Rev,digits)
if(min(Rev)>0) ans<-TRUE else(ans<-FALSE)
if(min(Im(v))>0) leigvec<-v else(leigvec<-Rev)
if(return.eigvec){
    return(list(ergodic=ans,eigvec=leigvec))
}
else{
    return(ans)
}}

