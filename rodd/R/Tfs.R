Tfs <- function(x, w, eta.f, eta.s, theta.f, theta.s)
{
    temp <- (eta.f(x, theta.f) - eta.s(x, theta.s))
    sum(w * temp * temp)
}
