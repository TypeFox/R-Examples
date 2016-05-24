`tyler.shape` <-
function(X,location=NULL,init=NULL,steps=Inf,eps=1e-6,maxiter=100,in.R=FALSE,print.it=FALSE,na.action=na.fail)
    {
    X<-na.action(X)
    if(!in.R) na.fail(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    
    if (is.numeric(location))
        {data.centered<-as.matrix(sweep(X,2,location,"-"))}
    else 
        {data.centered<-as.matrix(sweep(X,2,colMeans(X),"-"))}
        
    if(is.finite(steps)) maxiter<-Inf
    
    p<-dim(X)[2]
               
    if (p<2) stop("'X' must be at least bivariate")  
     
    center.ind<-apply(data.centered,1,setequal,y=rep(0,p))
    n.del<-sum(center.ind)
    
    if (n.del!= 0)
        {
        data.centered<-data.centered[center.ind==F,]
        if (n.del>1)
            {warning(paste(n.del ,"observations equal to the location center were removed"))}
        else
            {warning("One observations equal to the location center was removed")}
        }
    n<-dim(data.centered)[1]
    
    iter=0
    if (is.numeric(init)) V.0<-solve(init)
    else V.0<-solve(t(data.centered)%*%data.centered)
    differ=Inf
    
    if (!in.R)
    {
    while(TRUE)
        {
        if (any(iter>=steps,differ<eps)) break
        if (iter>=maxiter)
            {
             stop("maxiter reached without convergence")
            }
        # print(iter)
        
          sqrtV<-mat.sqrt(V.0)
          V.new<-sqrtV%*%(solve(sumsignout(data.centered%*%sqrtV)))%*%sqrtV
          V.new<-V.new/sum(diag(V.new))
          differ=frobenius.norm(V.new-V.0)
          V.0<-V.new
          iter=iter+1 
        }
    }
    else
    {
    while(TRUE)
        {
        if (any(iter>=steps,differ<eps)) break
        if (iter>=maxiter)
            {
             stop("maxiter reached without convergence")
            }
        # print(iter)
             V.new<-.tyler.step(V.0,data.centered,p,n)
            differ=frobenius.norm(V.new-V.0)
            V.0<-V.new
            iter=iter+1 
            }
    }
    if (print.it)
        {
        if (iter<steps) print(paste("convergence was reached after",iter, "iterations"))
        else print(paste("algorithm stopped after",steps, "iterations"))
        }
        
    V.shape<-solve(V.new)
    V<-V.shape/det(V.shape)^(1/p)
    
    colnames(V) <- colnames(X) 
    rownames(V) <- colnames(X) 
    return(V)
    
    }
