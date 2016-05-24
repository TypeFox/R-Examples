"mult" <-
function(n,r,p){
    exp( r*log(p)+(n-r)*log(1-p) 
        - lgamma(r+1) - lgamma(n-r+1) +lgamma(n+2) )
}

