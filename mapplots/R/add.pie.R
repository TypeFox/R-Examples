add.pie <-
function (z, x=0, y=0, labels=names(z), radius = 1, edges = 200, 
                     clockwise = TRUE, init.angle = 90, density = NULL, angle = 45,
                     col = NULL, border = NULL, lty = NULL, label.dist=1.1, ...) 
{
# modified pie() function from graphics library
# radius is scaled to units on the y axis
    if (!is.numeric(z) || any(is.na(z) | z < 0)) 
        stop("'z' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(z))
    else labels <- as.graphicsAnnot(labels)
    z <- c(0, cumsum(z)/sum(z))
    dz <- diff(z)
    nz <- length(dz)
    asp <- get.asp()
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nz)
    if (!is.null(border)) 
        border <- rep_len(border, nz)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nz)    
    angle <- rep(angle, nz)
    if (!is.null(density)) 
        density <- rep_len(density, nz)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = asp * radius * cos(t2p) +x , y = radius * sin(t2p) + y)
    }
    for (i in 1L:nz) {
        n <- max(2, floor(edges * dz[i]))
        P <- t2xy(seq.int(z[i], z[i + 1], length.out = n))
        polygon(c(P$x, 0+x), c(P$y, 0+y), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])    
        P <- t2xy(mean(z[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            #lines(c(P$x, (label.dist-0.05)*(P$x-x) + x),
            #      c(P$y, (label.dist-0.05)*(P$y-y)+y ) )
            text(label.dist*(P$x-x)+x , label.dist*(P$y-y)+y , labels[i], xpd = TRUE, 
                adj = ifelse(P$x-x < 0, 1, 0), ...)

        }
    }
}

