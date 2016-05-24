rmult.crm <-
function(clsize,lin.pred,cor.matrix,intercepts,link="probit")
{
 if(!is.numeric(clsize) | clsize < 2)
     stop("'clsize' must be greater than or equal to two")
 clsize <- as.integer(clsize) 
 lin.pred <- as.matrix(lin.pred)
 if(!is.numeric(lin.pred))
     stop("'lin.pred' must be a numeric")
 if(ncol(lin.pred)!=clsize) 
     stop("'lin.pred' must have ",clsize," columns")
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
 distr <- switch(link,"probit"="normal","logit"="logistic", "cloglog"="extreme","cauchit"="cauchit")
 if(!is.numeric(cor.matrix)) 
    stop("'cor.matrix' must be numeric")
 if(!is.matrix(cor.matrix))
    stop("'cor.matrix' must be matrix")
 if (!is.numeric(cor.matrix)) 
        stop("'cor.matrix' must be numeric")
    cor.matrix <- as.matrix(cor.matrix)
    dimcor <- clsize*(ncategories-1)
    if (ncol(cor.matrix) != dimcor) 
        stop("'cor.matrix' must be a ", dimcor, "x", dimcor, " matrix")
    if (!isSymmetric(cor.matrix)) 
        stop("'cor.matrix' must be symmetric")
    for (i in 1:clsize) {
        diag.index <- 1:(ncategories-1) + (i - 1) * (ncategories-1)
        cor.matrix[diag.index, diag.index] <- diag(1, ncategories-1)
    }
    if (any(cor.matrix > 1) | any(cor.matrix < -1)) 
        stop("all the elements of 'cor.matrix' must be on [-1,1]")
    if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0)) 
        stop("'cor.matrix' must be positive definite")
 err <- rnorta(R=R,cor.matrix=cor.matrix,distr=distr)
 lin.pred.extended <- t(apply(lin.pred,1,function(x) rep(x,each=ncategories-1)))
 U <- if(distr=="extreme") lin.pred.extended+err else -lin.pred.extended+err
 Ysim <- matrix(0,R,dimcor)
 for(i in 1:clsize) Ysim[,((i-1)*(ncategories-1)+1):(i*(ncategories-1))] <- cut(U[,((i-1)*(ncategories-1)+1):(i*(ncategories-1))],intercepts[i,],labels=FALSE)
 Ysim <- matrix(Ysim,ncol=ncategories-1,byrow=TRUE)
 for(i in 1:(ncategories-1)) Ysim[,i] <- ifelse(Ysim[,i]==i,i,ncategories)
 Ysim <- apply(Ysim,1,min)
 Ysim <- matrix(Ysim,R,clsize,byrow=TRUE)
 if(distr=="extreme") Ysim <- ncategories-Ysim+1
 list(Ysim=Ysim,correlation=cor.matrix,rlatent=err)
 }