NextPow2<-function (x) 
{
    a <- ceiling(log(x, 2))
    y <- 2^a
    y
}
