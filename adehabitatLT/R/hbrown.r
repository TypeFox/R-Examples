hbrown <- function(x)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class ltraj")

    so <- unlist(lapply(x, function(df) {
        df <- df[!is.na(df$x),]
        df <- df[!is.na(df$y),]
        dy <- df$y[-1]-df$y[-nrow(df)]
        dx <- df$x[-1]-df$x[-nrow(df)]
        dt <- unclass(df$date[-1])-unclass(df$date[-nrow(df)])
        df <- data.frame(x=df$x[-nrow(df)],y=df$y[-nrow(df)],dx=dx,dy=dy,dt=dt,
                     dist = sqrt(dx^2+dy^2))

        oo <- as.matrix(df[,c("dx","dy", "dt")])
        oo[,1] <- oo[,1]/sqrt(oo[,3])
        oo[,2] <- oo[,2]/sqrt(oo[,3])

        vc <- crossprod(scale(oo[,1:2],scale=FALSE))/nrow(oo)
        h <- sqrt(mean(diag(vc)))
        return(h)
    }))
    names(so) <- burst(x)
    return(so)
}
