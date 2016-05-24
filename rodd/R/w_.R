w_ <- function(alpha, w, k.max, k.min)
{
    w_ <- w
    w_[k.max] <- w_[k.max] + alpha
    w_[k.min] <- w_[k.min] - alpha 
    w_
}
