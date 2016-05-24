correlation.limits <-
function(n.BB,n.NN, prop.vec=NULL,coef.mat=NULL) {
 
    validation.bin(n.BB, prop.vec)
      
    if (missing(n.NN) == TRUE && !is.null(coef.mat)) {
        stop("Number of continuous variables is not specified !")
    } else
    if (n.NN > 0 && is.null(coef.mat))         {
        stop("Coefficient matrix is not specified while n.NN > 0 !")
    } else
    if (!is.null(coef.mat)) {
        if(n.NN == 0) {
        stop("Coefficient matrix is specified while n.NN=0")
        } else 
        if (n.NN > 0 && (ncol(coef.mat) != n.NN)) {
        stop("Dimension of coefficient matrix does not match the number of continuous variables! \n")
        }

     } #if


    if(!is.null(prop.vec)) {

    q.vec=(1-prop.vec)

    a=unlist(sapply(2:n.BB , function(i) sapply(1:(i-1), function(j) -sqrt((prop.vec[i]*prop.vec[j])/(q.vec[i]*q.vec[j])) )))
    b=unlist(sapply(2:n.BB , function(i) sapply(1:(i-1), function(j) -sqrt((q.vec[i]*q.vec[j])/(prop.vec[i]*prop.vec[j])) )))
    low.lim.b=apply(cbind(a,b),1,max)
    c=unlist(sapply(2:n.BB , function(i) sapply(1:(i-1), function(j) sqrt((prop.vec[i]*q.vec[j])/(q.vec[i]*prop.vec[j])) )))
    d=unlist(sapply(2:n.BB , function(i) sapply(1:(i-1), function(j) sqrt((q.vec[i]*prop.vec[j])/(prop.vec[i]*q.vec[j])) )))
    upp.lim.b=apply(cbind(c,d),1,min)

    sugcormat.b=diag(1,n.BB)
    sugcormat.b[lower.tri(sugcormat.b)]=low.lim.b
    sugcormat.b[upper.tri(sugcormat.b)]=upp.lim.b
    
    } #ifbinary


    if(!is.null(coef.mat)) {

    samples = 1e+05

    xmat=matrix(NA, nrow=samples, ncol=n.NN)
    for (i in 1:n.NN){
    x=as.vector(rnorm(samples))
    xx=cbind(1,x,x^2,x^3)
    xmat[,i]=xx%*%coef.mat[,i]
    }

    if(!is.null(prop.vec)){
    xmat=cbind(sapply(1:length(prop.vec),function(i) rbinom(samples,1,prop.vec[i])),xmat)
    }
        
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]   

    } #if

    if(!is.null(prop.vec) && is.null(coef.mat) ) {
    sugcormat= sugcormat.b
    diag(sugcormat)=NA
    } else 
    if( is.null(prop.vec) && !is.null(coef.mat) ) {
    sugcormat=diag(1,n.NN)
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    diag(sugcormat)=NA
    } else
    if( !is.null(prop.vec) && !is.null(coef.mat) ) {
    sugcormat=diag(1,(n.BB+n.NN))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[1:n.BB,1:n.BB]=sugcormat.b
    diag(sugcormat)=NA
    }
    
limitscor.mat=sugcormat

return(limitscor.mat)
}
