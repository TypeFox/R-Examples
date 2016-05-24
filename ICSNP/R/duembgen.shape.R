`duembgen.shape` <-
function(X,init=NULL,steps=Inf,eps=1e-6,maxiter=100,in.R=FALSE,na.action=na.fail,...)
    {
    X<-na.action(X)
    X<-as.matrix(X)
    
    n <- dim(X)[1] 
    p <- dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")  
    
    if (in.R==FALSE){
            if (is.null(init)) init <- cov(X)
            init <- shape.det(init)
            if(is.finite(steps)) maxiter <- Inf
            iter <- 0
            V <- init
            while(TRUE)
                {
                if(iter>=steps) break
                if(iter>=maxiter) stop("maxiter reached")
                iter <- iter + 1
                sqrtV <- mat.sqrt(V)
                V.new <- crossprod(sqrtV, SSCov(tcrossprod(X, solve(sqrtV)))) %*% sqrtV
                V.new <- shape.det(V.new)
                if(all(is.infinite(steps), base::norm(V.new-V, type="F")<eps)) break
                V<-V.new
                }
            duembgen <- V.new
            colnames(duembgen) <- colnames(X)
            rownames(duembgen) <- colnames(X)
            }
            
    
    if (in.R==TRUE){
        data2<-pair.diff(X)
        colnames(data2) <- colnames(X)
        duembgen<-tyler.shape(data2,0,init,steps,eps,maxiter,in.R,...)}
        
    return(duembgen)
    }


# Function to standardize a scatter matrix to have det 1, internal for duembgen.shape
shape.det <- function(M) M/det(M)^(1/dim(M)[2])

# Internal function for duembgen.shape
SSCov <- function(X)
{
d<-dim(X)
matrix(.C("sum_of_diff_sign_outers", as.double(X),as.integer(d), res=double(d[2]^2))$res,ncol=d[2],byrow=T)/(dim(X)[1]*(dim(X)[1]-1)/2)
}
