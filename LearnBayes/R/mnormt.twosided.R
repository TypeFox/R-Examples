mnormt.twosided <-
function (m0, prob, t, data) 
{
    xbar = data[1]
    n = data[2]
    h = data[3]
    num = 0.5 * log(n) - log(h) - 0.5 * n/h^2 * (xbar - m0)^2
    den = -0.5 * log(h^2/n + t^2) - 0.5/(h^2/n + t^2) * (xbar - 
        m0)^2
    bf = exp(num - den)
    post = prob * bf/(prob * bf + 1 - prob)
    return(list(bf = bf, post = post))
}

