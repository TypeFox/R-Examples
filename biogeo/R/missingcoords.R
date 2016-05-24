missingcoords <-
function (x, y) 
{
    xf <- getformat(x)
    xf1 <- matrix(unlist(str_locate(xf, "0")), ncol = 2, byrow = F)
    if (any(is.na(xf1[, 1]))) {
        f1 <- which(is.na(xf1[, 1]))
    }
    else {
        f1 <- {
        }
    }
    yf <- getformat(y)
    yf1 <- matrix(unlist(str_locate(yf, "0")), ncol = 2, byrow = F)
    if (any(is.na(yf1[, 1]))) {
        f2 <- which(is.na(yf1[, 1]))
    }
    else {
        f2 <- {
        }
    }
    f <- unique(c(f1, f2))
    return(f)
}
