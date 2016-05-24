#' Add a transient to a given mcmc sequence
#' @param X sequence to add the transient on
#' @param a last iteration of the constant transient part
#' @param b last iteratio of the transient
#' @param step transient step
#' @examples
#' require(codadiags)
#' x = AR1()
#' plot(x,type='l',col=rgb(.5,0,0,.5))
#' y = add.transient(x)
#' lines(y,col=rgb(0,0,0.5,.5))
#' transient.test(x)
#' transient.test(y)
add.transient <- function(X, a=100, b=a+100, step=-1) {
    for (i in 1:a){
        X[i] = X[i]+step
    }
    if (b > a+1)
        for (i in (a+1):b){
            X[i] = X[i] + (b-i)/(b-a) * step
        }
    return(X)
}
