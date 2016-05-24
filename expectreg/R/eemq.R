eemq <-
function (asy, ncp = 0, s = 1) 
{
    zz = 0 * asy
    for (k in 1:length(asy)) {
        root = function(z) peemq(z) - asy[k]
        z = uniroot(root, interval = c(-15, 15), tol = 1e-06)
        zz[k] = z$root * s + ncp
    }
    return(zz)
}
