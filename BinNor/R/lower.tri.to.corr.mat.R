lower.tri.to.corr.mat <-
function(corr.vec=NULL,d){
if(length(corr.vec)!=(d*(d-1)/2)){
stop("Vector of correlations is misspecified, dimension is wrong!\n")}
corr.mat = diag(d)
corr.mat [ lower.tri(corr.mat,diag=FALSE)]<-corr.vec
corr.mat = corr.mat + t(corr.mat)
diag(corr.mat)=1
return(corr.mat)
}

