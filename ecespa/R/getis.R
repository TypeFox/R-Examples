getis <-function (mippp, nx = 30, ny = 30, R = 10) 
{
   
    dataname <- deparse(substitute(mippp))
    cosagrid <- gridcenters(mippp$window, nx = nx, ny = ny)
    cosa1 <- ppp(x = mippp$x[1], y = mippp$y[1], window = mippp$window, 
        marks = "1")
    cosa2 <- setmarks(mippp, "2")
    cosa12 <- superimpose(cosa1, cosa2)
    cosa12$marks=factor(cosa12$marks)
    coso12 <- cosa12
    cosagrid <- as.data.frame(cosagrid)
    cosagrid <- rbind(cbind(x = mippp$x, y = mippp$y), cosagrid)
    klocalgrid <- apply(cosagrid, 1, function(x, cosa12 = coso12) {
        cosa12$x[1] <- x[[1]]
        cosa12$y[1] <- x[[2]]
        Kcross(unique(cosa12), i = "1", j = "2", r = c(0, R), 
            correction = "isotropic")$iso[2]
    })
    result <- list(x = cosagrid$x, y = cosagrid$y, klocal = klocalgrid[(1:mippp$n)], 
        klocalgrid = klocalgrid, R = R, nx = nx, ny = ny, dataname = dataname, 
        ppp = mippp)
    class(result) <- c("ecespa.getis", class(result))
    return(result)
}
