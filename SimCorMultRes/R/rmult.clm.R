rmult.clm <-
function(clsize,lin.pred,cor.matrix,intercepts,link="probit")
{
 if(!is.numeric(clsize) | clsize < 2)
     stop("'clsize' must be greater than or equal to two")
 clsize <- as.integer(clsize) 
 lin.pred <- as.matrix(lin.pred)
 if(!is.numeric(lin.pred))
     stop("'lin.pred' must be a numeric")
 if(ncol(lin.pred)!=clsize) 
     stop("the matrix 'lin.pred' must have ",clsize," columns")
 R <- nrow(lin.pred)
 if(!is.vector(intercepts) & !is.matrix(intercepts))
    stop("'intercepts' must be a vector or a matrix")
 if(!is.numeric(intercepts))
    stop("'intercepts' must be numeric")
 if(is.vector(intercepts)) {
   if(length(intercepts)==1)
     stop("'intercepts' must have at least 2 elements")
   if(any(diff(intercepts)<=0)) 
     stop("'intercepts' must be increasing") 
   ncategories <- length(intercepts)+1
   intercepts <- matrix(intercepts,clsize,ncategories-1,TRUE)
   intercepts <- cbind(-Inf,intercepts,Inf)
 } else {
   ncategories <- ncol(intercepts)+1
   intercepts <- cbind(-Inf,intercepts,Inf)
   for (i in 1:clsize){
     if(any(diff(intercepts[i,])<=0)) 
       stop("'intercepts' must be increasing at each row") 
   }
 }
 links <- c("probit","logit","cloglog","cauchit")
 if(!is.element(link,links)) 
   stop("'link' must be either 'probit','logit','cloglog' or'cauchit'") 
 distr <- switch(link,"probit"="normal","logit"="logistic",
                       "cloglog"="extreme","cauchit"="cauchit")
 if(!is.numeric(cor.matrix)) 
    stop("'cor.matrix' must be numeric")
 if(!is.matrix(cor.matrix))
    stop("'cor.matrix' must be matrix")
 if(ncol(cor.matrix)!=clsize | nrow(cor.matrix)!=clsize) 
    stop("'cor.matrix' must be a ",clsize,"x",clsize," matrix")
 if(!isSymmetric(cor.matrix)) 
    stop("'cor.matrix' must be a symmetric matrix") 
 if(any(diag(cor.matrix)!=1)) 
    stop("the diagonal elements of 'cor.matrix' must be one")
 if(any(cor.matrix>1) | any(cor.matrix< -1))
    stop("all the elements of 'cor.matrix' must be on [-1,1]")
 if(any(eigen(cor.matrix,symmetric=TRUE,only.values=TRUE)$values<=0))
    stop("'cor.matrix' must be positive definite")
 err <- rnorta(R=R,cor.matrix=cor.matrix,distr=distr)
 U <- if(distr=="extreme") lin.pred+err else -lin.pred+err
 Ysim <- matrix(0,R,clsize)
 for(i in 1:clsize) Ysim[,i] <- cut(U[,i],intercepts[i,],labels=FALSE)
 if(distr=="extreme") Ysim <- ncategories-Ysim+1
 list(Ysim=Ysim,correlation=cor.matrix,rlatent=err)
 }