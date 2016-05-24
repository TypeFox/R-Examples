cent.norm<-function (x,na.rm = FALSE) 
{
    x.cent <- sweep(x, 2, colMeans(x, na.rm = na.rm), "-")
    x.norm <- apply(x.cent^2, 2, sum, na.rm = na.rm )
    x.cent.norm <- sweep(x.cent, 2, sqrt(x.norm), "/")
    return(x.cent.norm)
}

