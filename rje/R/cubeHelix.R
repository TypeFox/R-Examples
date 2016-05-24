cubeHelix <-
function (n, start = 0.5, r = -1.5, hue = 1, gamma = 1) 
{
    M = matrix(c(-0.14861, -0.29227, 1.97294, 1.78277, -0.90649, 
        0), ncol = 2)
    lambda = seq(0, 1, length.out = n)
    l = rep(lambda^gamma, each = 3)
    phi = 2 * pi * (start/3 + r * lambda)
    t = rbind(cos(phi), sin(phi))
    out = l + hue * l * (1 - l)/2 * (M %*% t)
    out = pmin(pmax(out, 0), 1)
    out = apply(out, 2, function(x) rgb(x[1], x[2], x[3]))
    return(out)
}
