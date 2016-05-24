complexity <- function(method, exchangeable, r, p, display = TRUE) {
    if (method == "Bonferroni") {
        if (exchangeable) {
            true.complexity <- p - r + 1
        } else {
        # I have a problem installing package hypergeo
        # cat(paste("Complexity: ", choose(p, r) * hypergeo(1, r - p, r + 1, -1), "\n"))
            if (p - r + 1 > r) {
                true.complexity <- sum(choose(p, (p - r + 1):p))
            } else {
                true.complexity <- sum(choose(p, r:p))
            }            
        }
    }
    if (method == "Hochberg") {
        if (exchangeable) {
            true.complexity <- choose(2 * p - r + 1, p - r + 1) - choose(2 * p - r + 1, p - r)
        } else {
            if (display) cat(paste("Order of complexity: ", (p - r + 2) ^ p, "\n"))
            w <- 1:(p - r + 1)
            true.complexity <- sum(apply(.ind.sumastar(w, p)$Mat, FUN = function(a.vec, p) .coeff.mnom(diff(c(0, a.vec, p))), MARGIN = 2, p = p))
        }
    }
    if (method == "Holm") {
        if (exchangeable) {
            true.complexity <- choose(p + r, r) - choose(p + r, r - 1)
        } else {
            if (display) cat(paste("Order of complexity: ", (r + 1) ^ p, "\n"))
            w <- 1:r
            true.complexity <- sum(apply(.ind.sumastar(w, p)$Mat, FUN = function(a.vec, p) .coeff.mnom(diff(c(0, a.vec, p))), MARGIN = 2, p = p))
        }
    }
    if (display) cat(paste("True complexity: ",true.complexity , "\n"))
    return(true.complexity)
}



    
        
