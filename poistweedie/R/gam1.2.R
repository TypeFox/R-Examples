gam1.2 <-
function (y, lambda) 
{
    gamma <- double(1)
    gamma <- exp((y + lambda - 1/2) * log(y + lambda) + (1/(12 * 
        (y + lambda))) - (y + 1/2) * log(y) - (lambda - 1/2) * 
        log(lambda) - (1/2) * log(2 * pi) - 1/(12 * lambda))
    gamma
}

