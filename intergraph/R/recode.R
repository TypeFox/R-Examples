# given a vector and a recode matrix replace the values in 'x' that match in
# 'm[,1]' with the corresponding replacements from 'm[,2]'

recode <- function(x, mat)
{
    i <- match(x, mat[,1])
    i2 <- which(!is.na(i))
    i <- i[i2]
    rval <- x
    rval[i2] <- mat[,2][i]
    rval
}
