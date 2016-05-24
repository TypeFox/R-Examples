ellipse = function(xcen = 0, ycen = 0, a = 10, b = 5, e = 1-b/a, pa = 0){
    # calculates ellipse x,y coordinates for plotting
    # Xcenter, Ycenter, semi-major axis, ellipticity(1-b/a), position angle (x+ ACW) (degrees)
    alpha = seq(0, 2 * pi, length=1000)
    pa = (pi/180)*pa
    b = (1-e)*a
    x = xcen + a * cos(alpha) * cos(pa) - b * sin(alpha) * sin(pa)
    y = ycen + a * cos(alpha) * sin(pa) + b * sin(alpha) * cos(pa)
    return(list(x=x,y=y))
}

