ARI <- function (x=NULL, y=NULL)
{
    x <- as.vector(x)
    y <- as.vector(y)
    
    tab <- table(x,y)
    
    ARI= classAgreement(tab)$crand
    return(ARI)
}
