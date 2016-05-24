work.corr <-
function(object, digits = 3){
 if(is.null(object$poly$poly) == TRUE){
  xsmat <- as.matrix(smat(coeff = object$coefficients[1:(object$categories-1)])$smat)
 } else {
  xsmat <- as.matrix(smat(coeff = object$poly$polycuts$coeff)$smat)
 }
 if(object$corr.mod != "independence"){
  xicmat <- cmat(ctimes = object$times, alpha = object$alpha, corrmod = object$corr.mod, 
        diffmeth = object$diffmeth, h = object$fit.opt[5])$icmat
  xcmat <- as.matrix(solve(xicmat)) } else {
  xcmat <- as.matrix(Matrix::Diagonal(length(object$times)))
 }
 wcor <- kronecker(xcmat, xsmat)
 rownames(wcor) <- colnames(wcor) <- 1:object$max.id
 return(round(wcor, digits = digits))
}
