"ctsub"<-
function(x, y, z)
{
    ##
    ## for function defined by (x(1),y(1))...(x(n),y(n)), compute
    ## integral from 0 to z(i) and put it in ans(i), for i=1,2,..n
    ## uses the trapezoid rule
    ##
    ## used by boott
    
    n <- length(z)
    ans <- rep(0,n)
    for(i in 1:n)
    {
        if(z[i]<= x[1]) {ans[i] <- (z[i]-x[1])*y[1]}
        else {
            j <- 1
            ans[i] <- 0
            while((j<=n) & (z[i]>x[j]) ){
                if(j > 1){
                    ans[i] <- ans[i]+(x[j]-x[j-1])*(y[j]+y[j-1])/2
                }
                j <- j+1
            }
            if(z[i] <= x[n]){
                ans[i] <-
                    ans[i]+.5*(z[i]-x[j-1])*(2*y[j-1]+(z[i]-x[j-1])*(y[j]-y[j-1])/(x[j]-x[j-1])) 
            }
            else { ans[i] <- ans[i]+(z[i]-x[n])*y[n] }
        }
    }
    
    return(ans)
}
