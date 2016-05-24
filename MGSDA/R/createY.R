.createY <-function(ytrain){
    n=length(ytrain)
    G=max(ytrain)
    # Create matrix Z
    Z=matrix(0,n,G)
    for (g in 1:G){
        Z[ytrain==g,g]=1
    }
    # Create cumulative sums
    cumsum=rep(sum(ytrain==1),G)
    for (i in 2:G){
        cumsum[i]=cumsum[i-1]+sum(ytrain==i)
    }
    
    # Create matrix H
    H=matrix(0,G,G-1)
    for (i in 1:(G-1)){
        # initialize each column
        H[1:i,i]=sqrt(sum(ytrain==i+1)/(cumsum[i]*cumsum[i+1]))
        H[i+1,i]=-sqrt(cumsum[i]/(cumsum[i+1]*sum(ytrain==i+1)))
    }
    
    # Create Y
    return(sqrt(n)*(Z%*%H))
}