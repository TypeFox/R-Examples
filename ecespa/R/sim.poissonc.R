`sim.poissonc` <-
function (x.ppp, rho, sigma) 
{
    mirpoispp = function(lambda, win) {
        n = lambda * area.owin(win)
        x = runif(n, win$xrange[1], win$xrange[2])
        y = runif(n, win$yrange[1], win$yrange[2])
        return(ppp(x, y, window = win))
    }
    ventana = x.ppp$window
    npadres = rho * area.owin(ventana)
    media = x.ppp$n/npadres
    padres = mirpoispp(lambda = rho, win = ventana)
    padres.pois = rpois(npadres, media)
    hijosx = rep(padres$x, padres.pois)
    hijosy = rep(padres$y, padres.pois)
    desviacionesx = rnorm(sum(padres.pois), mean = 0, sigma)
    desviacionesy = rnorm(sum(padres.pois), mean = 0, sigma)
    hijosx = hijosx + desviacionesx
    hijosy = hijosy + desviacionesy
    for (i in 1:length(hijosx)) {
        if (hijosx[i] < ventana$x[1]) 
            hijosx[i] = ventana$x[2] - (abs(hijosx[i]) - abs(ventana$x[1]))
        if (hijosx[i] > ventana$x[2]) 
            hijosx[i] = ventana$x[1] + (abs(hijosx[i]) - abs(ventana$x[2]))
        if (hijosy[i] < ventana$y[1]) 
            hijosy[i] = ventana$y[2] - (abs(hijosy[i]) - abs(ventana$y[1]))
        if (hijosy[i] > ventana$y[2]) 
            hijosy[i] = ventana$y[1] + (abs(hijosy[i]) - abs(ventana$y[2]))
    }
    return(ppp(x = hijosx, y = hijosy, window = ventana))
}

