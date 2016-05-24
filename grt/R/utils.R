old2new_par <- function(x)
{
    if(!inherits(x, "glcStruct")) 
        stop("x not of class 'glcStruct'")
    if(is.null(x$coeffs) & !is.null(x$angle))
       angle <- x$angle
    else 
       angle <- cart2angle(x$coeffs)
    new_params <- structure(
       list(noise = x$noise,
       angle = angle, 
       bias = x$bias),
       class=c("glcStruct", "list"))
    new_params
}

new2old_par <- function(x)
{
    if(!inherits(x, "glcStruct")) 
        stop("x not of class 'glcStruct'")
    if(!is.null(x$coeffs) & is.null(x$angle))
       coeffs <- x$coeffs
    else 
       coeffs <- angle2cart(x$angle)
    old_params <- structure(
       list(noise = x$noise,
       coeffs = coeffs, 
       bias = x$bias),
       class=c("glcStruct", "list"))
    old_params
}

cart2angle <- function(cart)
{
    x <- as.vector(cart)
    if(length(x) == 1) stop("no angle for 1D")
    else if(length(x) == 2)
    {
        x <- x / sqrt(sum(x^2))
	angles <- acos(x[1])
	if(x[2]<0) angles <- 2*pi - angles
    }
    else if(length(x) == 3)
    {
        x <- x / sqrt(sum(x^2))
	angles <- c(asin(x[3]), acos(x[1]/cos(asin(x[3]))) )
	if(x[2]<0) angles[2] <- 2*pi - angles[2]
    }
    else stop("length(x) should be between 1 and 3")
    angles
}

angle2cart <- function(angle)
{
    x <- as.vector(angle)
    if(length(x) == 1) cart <- c(cos(x), sin(x))
    else if(length(x) == 2)
        cart <- c(cos(x[1])*cos(x[2]), cos(x[1])*sin(x[2]), sin(x[1]))
    else stop("length(x) should be between 1 and 2")
    cart
}

unscale <- function(x)
{
    if(!is.null(attr(x,"scaled:scale")))
        x <- sweep(x, FUN="*", 2, attr(x,"scaled:scale"))
    if(!is.null(attr(x,"scaled:scale")))
        x <- sweep(x, FUN="+", 2, attr(x,"scaled:center"))
    attr(x,"scaled:scale") <- NULL
    attr(x,"scaled:center") <- NULL
    x
}