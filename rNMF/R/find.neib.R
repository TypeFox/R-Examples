## Given a line number in a vectorized column i, find neighboring points of it in a non-vectorized matrix.
## p1 = the number of rows of the individual image
## n1 = the number of columns of the individual image
## zeta = the matrix of points to keep
## A = the original image
## Not used so far.
find.neib = function(i, p1, n1, zeta, A, k){
    step = 1
    x = find.x(i,p1)
    y = find.y(i,p1)
    while(TRUE){
        p.range = (i - min(step, x - 1)) : (i + min(step, p1 - x)) 
        n.range = (-min(step, y - 1)) : (min(step, n1 - y))
        to.check = c(rep(p.range, length(n.range)) + rep(n.range, each = length(p.range)) * p1)
        if(all(!zeta[to.check,])){
            step = step + 1
        }else{
            return(rep(mean(A[to.check,][zeta[to.check,]]), k))
        }
    }
}
