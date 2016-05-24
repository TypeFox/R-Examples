gam1.1 <-
function (y, lambda) 
{
    gamma <- double(1)
    if (y == 0) {
        gamma <- 1
    }
    else {
        gamma <- 1
        for (i in 0:(y - 1)) {
            gamma <- (gamma * (lambda + i))/(i + 1)
        }
    }
    gamma
}

