## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

projected_gradient_norm <- function(l,u,x,g){
    norm=0
    for(i in 1:length(x)){
        if(g[i]<0){
            gi=max(x[i]-u[i],g[i])
        } else {
            gi=min(x[i]-l[i],g[i])
        }
        norm=max(norm,abs(gi))
    }
    return( norm)
}
