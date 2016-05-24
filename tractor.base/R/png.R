writePng <- function (data, fileName, aspectRatio = 1)
{
    if (length(dim(data)) == 3)
    {
        data <- data[,dim(data)[2]:1,]
        data <- aperm(data, c(2,1,3))
    }
    else if (length(dim(data)) == 2)
        data <- t(data[,dim(data)[2]:1])
    
    data[data<0] <- 0
    data[data>1] <- 1
    
    png::writePNG(data, ensureFileSuffix(fileName,"png"), asp=aspectRatio)
}
