plotcdf2 <- function (x, y, f, xaxe, yaxe, col = NULL, border = FALSE, Nxy = 200, 
    theme = "0") 
{
    if (length(f) > 1) {
        xi = sort(x)
        yj = sort(y)
        k = length(x)
        l = length(y)
    }
    else {
        xi = as.numeric(levels(as.factor(x)))
        yj = as.numeric(levels(as.factor(y)))
        f = table(x, y)
        k = length(xi)
        l = length(yj)
    }
    if (sum(sum(f)) > 1) {
        f = f/sum(sum(f))
    }
    F = matrix(0,ncol=l,nrow=k)
    F[1, ] = cumsum(f[1, ])
    F[, 1] = cumsum(f[, 1])
    for (i in 2:k) {
        for (j in 2:l) {
            F[i, j] = f[i, j] + F[i - 1, j] + F[i, j - 1] - F[i - 
                1, j - 1]
        }
    }
    deltax = (max(xi) - min(xi))/Nxy
    deltay = (max(yj) - min(yj))/Nxy
    x = seq(min(xi) - deltax, max(xi) + deltax, deltax)
    y = seq(min(yj) - deltay, max(yj) + deltay, deltay)
    n1 = length(x)
    n2 = length(y)
    z = matrix(rep(0, n1 * n2), ncol = n2)
    for (i in 1:n1) {
        for (j in 1:n2) {
            i1 = (x[i] >= xi)
            i2 = (y[j] >= yj)
            if (sum(i1) == 0 | sum(i2) == 0) {
                z[i, j] = 0
            }
            if (sum(i1) >= k & sum(i2) >= l) {
                z[i, j] = 1
            }
            if (sum(i1) >= k & sum(i2) < l & sum(i2) > 0) {
                z[i, j] = F[k, sum(i2)]
            }
            if (sum(i1) < k & sum(i2) >= l & sum(i1) > 0) {
                z[i, j] = F[sum(i1), l]
            }
            if (sum(i1) < k & sum(i2) < l & sum(i1) > 0 & sum(i2) > 
                0) {
                z[i, j] = F[sum(i1), sum(i2)]
            }
        }
    }
    if (is.null(col)) {
        nrz <- nrow(z)
        ncz <- ncol(z)
        jet.colors <- colorRampPalette(c("blue", "red"))
        if (theme == "1") {
            jet.colors <- colorRampPalette(c("#BDFF00", "#FF00BD", 
                "#00BDFF"))
        }
        if (theme == "2") {
            jet.colors <- colorRampPalette(c("#FF8400", "#8400FF", 
                "#00FF84"))
        }
        if (theme == "3") {
            jet.colors <- colorRampPalette(c("#84FF00", "#FF0084", 
                "#0084FF"))
        }
        if (theme == "bw") {
            jet.colors <- function(nbcols) {
                gray(seq(0.1, 0.9, length.out = nbcols))
            }
        }
        nbcol <- 100
        color <- jet.colors(nbcol)
        zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, 
            -ncz]
        facetcol <- cut(zfacet, nbcol)
        persp(x, y, z, theta = -30, phi = 15, col = color[facetcol], 
            shade = 0.15, main = "St\u00E9r\u00E9ogramme des deux variables", 
            xlab = xaxe, ylab = yaxe, zlab = "", cex.axis = 0.75, 
            ticktype = "detailed", border = border)
    }
    else {
        persp(x, y, z, theta = -30, phi = 15, col = col, shade = 0.15, 
            main = "St\u00E9r\u00E9ogramme des deux variables", xlab = xaxe, 
            ylab = yaxe, zlab = "", cex.axis = 0.75, ticktype = "detailed", 
            border = border)
    }
    invisible(list(F=F,z=z,x=x,y=y))
}
