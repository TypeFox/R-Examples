quadtrafo <- function(e, f=NULL, g=NULL, h=NULL)
{
    ## assumed to require() for all other quad**** function
    require(scatterplot3d)
    if (is.matrix(e)) dat <- e
    else if (is.null(f)) dat <- t(e)
        else dat <- cbind(e,f,g,h)
    result <- dat %*%   
        matrix(c(1, 0.5, 0.5, 0, 
                0, sqrt(3)/2, sqrt(3)/6, 0, 
                0, 0, sqrt(6)/3, 0), nrow = 4)
                
    colnames(result) <- c("x", "y", "z")
    return(result)
}


quadlines <- function(e, f=NULL, g=NULL, h=NULL, sp, ...)
{
  result <- quadtrafo(e,f,g,h) # requires scatterplot3d
  sp$points3d(result[,1], result[,2], result[,3], type="l", ...)
  invisible(result)
}


quadpoints <- function(e, f=NULL, g=NULL, h=NULL, sp, ...)
{
  result <- quadtrafo(e,f,g,h) # requires scatterplot3d
  sp$points3d(result[,1], result[,2], result[,3], type="p", ...)
  invisible(result)
}


quadplot <- function(e=NULL, f=NULL, g=NULL, h=NULL, angle=75, scale.y=0.6, 
                    label=1:4, labelcol=rainbow(4), labelpch=19, 
                    labelcex=1.5, main="", 
                    s3d.control = list(), simplex.control = list(), 
                    legend.control = list(), ...)
{

    corners <- quadtrafo(diag(4)) # requires scatterplot3d
    s3d <- do.call("scatterplot3d", 
        c(list(0.5, 0.2886751, 0.2041241, type="n", 
            xlim=range(corners[,1]), ylim=range(corners[,2]), zlim=range(corners[,3]),
            axis = FALSE, grid=FALSE, angle = angle, scale.y = scale.y, main = main),
        s3d.control)
    )

    ## we don't want to waste too much space, 
    ## hence setting optimized "usr" values into the environment of s3d:
    usr <- as.vector(sapply(s3d$xyz.convert(corners[-2,]), range))
    opar <- par(usr = usr, xpd = NA)
    assign("usr", usr, envir = environment(s3d$points3d))
    on.exit(par(opar))

    ## outer simplex:
    do.call("quadlines", 
        c(list(e = diag(4)[c(1:4,1,3,2,4),], sp = s3d), 
        simplex.control)
    )
    do.call("quadpoints", 
        c(list(e = diag(4), sp = s3d, pch = labelpch, col = labelcol, cex = labelcex)))
    ## legend:
    do.call("legend", 
        c(list(usr[1], usr[4], legend = label, col = labelcol, pch = labelpch, cex = labelcex),
        legend.control)
    )
    ## points:
    if (!is.null(e)) 
        quadpoints(e, f, g, h, sp = s3d, ...)
    return(s3d)
}
