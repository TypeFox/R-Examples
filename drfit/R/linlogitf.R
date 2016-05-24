linlogitf <- function(x,k,f,mu,b)
{
    k*(1 + f*x) / (1 + ((2*f*(10^mu) + 1) * ((x/(10^mu))^b)))
}
