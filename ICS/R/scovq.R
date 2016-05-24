scovq <- function(x,y,q1=0,q2=0.5, pos=TRUE, type=7, method = "unbiased", na.action=na.fail, check=TRUE)
    {
    if (check){
    nr <- nrow(x)
    if (nr != length(y)) stop("'x' and 'y' have different lengths")
    if(q1<0 | q1>=q2 | q2>1) stop("'q1' and 'q2' must be between 0 and 1 and  'q1' must be smaller than 'q2'")
    X <- data.frame(y)
    X$x <-as.matrix(x)
    X <- na.action(X)
    y <- X$y
    x <- X$x
    if(!all(sapply(x, is.numeric))) stop("'x' must be numeric")
    if(!is.numeric(y)) stop("'y' must be numeric")
    }
    x <- as.matrix(x)
    y.quant <- quantile(y,probs=c(q1,q2), type=type)
    wt.y <- ifelse(y < y.quant[1] | y > y.quant[2], 1, 0)
    if(!pos) wt.y <- 1 - wt.y
    wt.y <- wt.y/mean(wt.y)
    cov.wt(x, wt.y, method=method)$cov
    } 
 
