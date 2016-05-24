`runifsphere` <-
function(n,p)
    {
    p<-as.integer(p)
    if(!is.integer(p)) stop("p must be an integer larger or equal than 2")
    if(p<2) stop("p must be an integer larger or equal than 2")
    
    Mnormal <- matrix(rnorm(n*p,0,1),nrow=n)
    rownorms <- sqrt(rowSums(Mnormal^2))
    
    unifsphere <- sweep(Mnormal,1,rownorms, "/")
    return(unifsphere)
    }

