correlation.limits <-
function(n.P, n.B, n.O, lambda.vec=NULL, prop.vec=NULL, prop.list=NULL) {
 
    validation.bin(n.B, prop.vec)
    validation.ord(n.O, prop.list)
  
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

    
    if(!is.null(prop.list)) {
    samples = 1e+05
    xmat3=sapply(1:length(prop.list), function(i) apply(as.matrix(runif(samples),samples,1), 1, function(x) length(which(prop.list[[i]]<x))+1))
   
    sxmat=apply(xmat3,2,sort)
    upp.lim.o=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim.o=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]   
    
    sugcormat.o=diag(1,n.O)
    sugcormat.o[lower.tri(sugcormat.o)]=low.lim.o
    sugcormat.o[upper.tri(sugcormat.o)]=upp.lim.o
   
    } #ifordinal

    if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(prop.list) ) {
    sugcormat=sugcormat.p
    diag(sugcormat)=NA
    } else 
    if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list) ) {
    sugcormat=sugcormat.b
    diag(sugcormat)=NA
    } else 
    if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list) ) {
    sugcormat=sugcormat.o
    diag(sugcormat)=NA
    } else
    if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list)) {
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
    if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list)) {
    xmat=cbind(xmat1,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.P+n.O))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    diag(sugcormat)=NA
    } else
    if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list)) {
    xmat=cbind(xmat2,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.B+n.O))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[1:n.B,1:n.B]=sugcormat.b
    diag(sugcormat)=NA
    } else
    if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list)) {
    xmat=cbind(xmat1,xmat2,xmat3)
    sxmat=apply(xmat,2,sort)
    upp.lim=cor(sxmat)[col(cor(sxmat)) > row(cor(sxmat))] 
    rsxmat=apply(sxmat,2,rev)
    low.lim=cor(sxmat,rsxmat)[col(cor(sxmat,rsxmat)) < row(cor(sxmat,rsxmat))]  
    sugcormat=diag(1,(n.P+n.B+n.O))
    sugcormat[lower.tri(sugcormat)]=low.lim
    sugcormat[upper.tri(sugcormat)]=upp.lim
    sugcormat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=sugcormat.b
    #sugcormat[(n.P+n.B+1):(n.P+n.B+n.O),(n.P+n.B+1):(n.P+n.B+n.O)]=sugcormat.o
    diag(sugcormat)=NA
    } 

limits.corr.mat=sugcormat

return(limits.corr.mat)
}
