"factab" <- 
function (object) 
{
    if(class(object) != "car")
    stop("object must be car\n")
    call <- match.call()
    par <- object$phi
    p <- object$order
    k <- object$scale
    w <- sort(round(polyroot(rev(c(1, par[1:p]))),4))
    r <- -k * (1 - w)/(1 + w)
    sortind <- sort.list(abs(Re(r)), decreasing = TRUE)
    r <- rev(r[sortind]) #changed  
    car.freq <- abs(Im(r))/(2 * pi)
    structure(list(call=call,root=r,freq=car.freq),
              class="factab")
}

print.factab <- 

function (x, digits = 3, ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nCharacteristic root of original parameterization in alpha", "\n\n")
    root <- drop(round(x$root, digits=digits))
    freq <- drop(round(x$freq, digits=digits))
    names(root) <- seq(length = length(x$freq))
    print.default(root, print.gap = 2)
    cat("\nFrequency", "\n\n")
    names(freq) <- seq(length = length(x$freq))
    print.default(freq, print.gap = 2)
    invisible(x)
  }
