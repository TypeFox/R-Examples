plotsurf <- function(data, x=2, y=3, z=4, mode=1, ticktype="detailed", expand=0.75, d=100, theta=-30, phi=25, ...)
{
    if(!is.data.frame(data))
        data <- read.table(data, header=TRUE)

    if(is.character(x))
        x <- match(x, names(data))
    if(is.character(y))
        y <- match(y, names(data))
    if(is.character(z))
        z <- match(z, names(data))

    x <- data[, x]
    y <- data[, y]
    z <- data[, z]

    ## interp is from package akima, but this package does not have a namespace,
    ## so we have to load the whole package (see DESCRIPTION).
    if(!isTRUE(akimaPermitStatus()))
       stop("The licence for akima is not Free or Open Source and explicitly forbids commercial use.\nIf you want to use akima anyway, call akimaPermit() (see help(akimaPermit) for details.")
    data <- akima::interp(x, y, z)
    x <- data$x
    y <- data$y
    z <- data$z

    if(mode==1)
        persp(x, y, z, ticktype=ticktype, expand=expand, d=d, theta=theta, phi=phi, ...)
    else if(mode==2)
        image(x, y, z, ...)
    else if(mode==3)
        contour(x, y, z, ...)

    return(invisible())
}

akimaPermit <- function ()
{
    if ("akima" %in% .packages(all.available = TRUE))
        assign(".akimaStatus", TRUE, envir=.BAYESX_CACHE)
    akimaPermitStatus()
}

akimaPermitStatus <- function ()
  {
  get(".akimaStatus", envir=.BAYESX_CACHE)
  }
