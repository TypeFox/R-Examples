diagmod <-
function (x) 
{
    dimx <- length(x)
    ans <- matrix(0, dimx, dimx)
    ans[1 + 0:(dimx - 1) * (dimx + 1)] <- x
    ans
}

