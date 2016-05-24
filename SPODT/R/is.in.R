is.in <- function(x1,y1,x2,y2,a1,b1){
        x.min=min(x1,x2); x.max=max(x1,x2)
        y.min=min(y1,y2); y.max=max(y1,y2)
        if(x1 != x2){
            b <- (y2 - y1)*a1/(x2 -x1) + (y1 - (y2 - y1)*x1/(x2 - x1))
            in.seg <- is.logical(all.equal(b,b1))
        }
        else{
            in.seg <- is.logical(all.equal(a1,x1))
        }
        res.log <- (in.seg)&(x.min<a1)&(a1<x.max)&(y.min<b1)&(b1<y.max)
        return(res.log)
}
