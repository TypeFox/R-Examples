correlation.limits <-
function(n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat=NULL) {
 
    validation.bin(n.B, prop.vec)
  
    if (missing(n.P) == TRUE && !is.null(lambda.vec)) {
        stop("Number of Poisson variables is not specified!")
    } else
    if (n.P > 0 && is.null(lambda.vec))         {
        stop("Lambda vector is not specified while n.P > 0!")
    } else
    if (!is.null(lambda.vec)) {
        if(n.P == 0) {
        stop("Lambda vector is specified while n.P=0!")
        } else 
        if (n.P > 0 && (length(lambda.vec) != n.P)) {
        stop("Length of lambda vector does not match the number of Poisson variables! \n")
        } else
        errorCount1=0
        for (i in 1:length(lambda.vec)){
        if(lambda.vec[i] <= 0) {
        cat("\n Lambda for Poisson variable",i,"must be greater than '0'!","\n")
        errorCount1 = errorCount1 + 1
        cat("\n")
        } #if
        } #for 
        if (errorCount1 > 0) {
        stop("Range violation occurred in the lambda vector!")
        }#if
    } #if

    if (missing(n.C) == TRUE && !is.null(coef.mat)) {
        stop("Number of continuous variables is not specified!")
    } else
    if (n.C > 0 && is.null(coef.mat))         {
        stop("Coefficient matrix is not specified while n.C> 0!")
    } else
    if (!is.null(coef.mat)) {
        if(n.C == 0) {
        stop("Coefficient matrix is specified while n.C=0!")
        } else 
        if (n.C > 0 && (ncol(coef.mat) != n.C)) {
        stop("Dimension of coefficient matrix does not match the number of continuous variables! \n")
        }

     } #if


    if(!is.null(lambda.vec)) {
   
    samples=1e+05
    xmat1=sapply(1:length(lambda.vec),function(i) rpois(samples,lambda.vec[i]))

    sxmat=apply(xmat1,2,sort)
    upp.lim.p=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim.p=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]   

    sugcormat.p=diag(1,n.P)
    sugcormat.p[lower.tri(sugcormat.p)]=low.lim.p
    sugcormat.p[upper.tri(sugcormat.p)]=upp.lim.p
      
    } #ifpoisson

    if(!is.null(prop.vec)) {

    q.vec=(1-prop.vec)

    a=unlist(sapply(2:n.B , function(i) sapply(1:(i-1), function(j) -sqrt((prop.vec[i]*prop.vec[j])/(q.vec[i]*q.vec[j])) )))
    b=unlist(sapply(2:n.B , function(i) sapply(1:(i-1), function(j) -sqrt((q.vec[i]*q.vec[j])/(prop.vec[i]*prop.vec[j])) )))
    low.lim.b=apply(cbind(a,b),1,max)
    c=unlist(sapply(2:n.B , function(i) sapply(1:(i-1), function(j) sqrt((prop.vec[i]*q.vec[j])/(q.vec[i]*prop.vec[j])) )))
    d=unlist(sapply(2:n.B , function(i) sapply(1:(i-1), function(j) sqrt((q.vec[i]*prop.vec[j])/(prop.vec[i]*q.vec[j])) )))
    upp.lim.b=apply(cbind(c,d),1,min)

    samples = 1e+05
    xmat2=sapply(1:length(prop.vec),function(i) rbinom(samples,1,prop.vec[i]))

    sugcormat.b=diag(1,n.B)
    sugcormat.b[lower.tri(sugcormat.b)]=low.lim.b
    sugcormat.b[upper.tri(sugcormat.b)]=upp.lim.b
      
    } #ifbinary

    
    if(!is.null(coef.mat)) {

    samples = 1e+05

    xmat3=matrix(NA, nrow=samples, ncol=n.C)
    for (i in 1:n.C){
    x=as.vector(rnorm(samples))
    xx=cbind(1,x,x^2,x^3)
    xmat3[,i]=xx%*%coef.mat[,i]
    }

    sxmat=apply(xmat3,2,sort)
    upp.lim.n=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim.n=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]   

    sugcormat.n=diag(1,n.C)
    sugcormat.n[lower.tri(sugcormat.n)]=low.lim.n
    sugcormat.n[upper.tri(sugcormat.n)]=upp.lim.n
    }

    if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(coef.mat) ) {
    sugcormat=sugcormat.p
    diag(sugcormat)=NA
    } else 
    if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat) ) {
    sugcormat=sugcormat.b
    diag(sugcormat)=NA
    } else 
    if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat) ) {
    sugcormat=sugcormat.n
    diag(sugcormat)=NA
    } else
    if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat)) {
    xmat=cbind(xmat1,xmat2)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.P+n.B))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=sugcormat.b
    diag(sugcormat)=NA
    } else
    if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat)) {
    xmat=cbind(xmat1,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.P+n.C))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    diag(sugcormat)=NA
    } else
    if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat)) {
    xmat=cbind(xmat2,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.B+n.C))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[1:n.B,1:n.B]=sugcormat.b
    diag(sugcormat)=NA
    } else
    if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat)) {
    xmat=cbind(xmat1,xmat2,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.P+n.B+n.C))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=sugcormat.b
    diag(sugcormat)=NA
    } 

limits.corr.mat=sugcormat

return(limits.corr.mat)
}
