heatmap <- function(ldlasso.obj, s2.indx = NA){

#par(mfrow = c(2,1) )

M <- abs(ldlasso.obj$cp.obj$beta0.mat)[ldlasso.obj$cp.obj$s2.vec <= 1e6,]
N1 <- nrow(M); N2 <- ncol(M);

scale = 1e1
image(1:ncol(M), 1:nrow(M), t(M), col=grey((1:scale)/(1+scale))[scale:1], breaks=c(seq( 1e-6, max(M), length.out = 1+scale )), zlim = c(1e-6,max(M)), xlab="", ylab="", axes=FALSE, cex=1)
if(!is.na(s2.indx)){
  abline( h = (1:length(ldlasso.obj$cp.obj$s2.vec))[s2.indx], lty = 3 )
}
abline( h = ldlasso.obj$s2star, lty = 2 )
axis(1)
box()

#fines = 1e3
#M <- ifelse( M > 1e-6, 1, 0 )
#image(1:ncol(M), 1:nrow(M), t(M), col=grey(c(1,0)), breaks=c(0,1/2,1) , xlab="SNP", ylab="", axes=FALSE, cex=1)

#if(!is.na(s2.indx)){
#  abline( h = (1:length(ldlasso.obj$cp.obj$s2.vec))[s2.indx], lty = 3 )
#}
#abline( h = ldlasso.obj$s2star, lty = 2 )
#axis(1)
#box()

}
