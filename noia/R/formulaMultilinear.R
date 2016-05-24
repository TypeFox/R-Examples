formulaMultilinear <-
function (nloc = 2, max.level = 2, max.dom = 2, e.unique = FALSE) 
{
    a <- noia::effectsNames[2]
    d <- noia::effectsNames[3]
    f <- paste("phen ~ X[[\"", effNames(nloc = nloc), "\"]]*R", 
        sep = "")
    for (l1 in 1:nloc) {
        add <- effNames(c(a), c(l1), nloc)
        dom <- effNames(c(d), c(l1), nloc)
        f <- paste(f, " + ", "X[[\"", add, "\"]]*a", l1, sep = "")
        if (max.dom > 0) {
            f <- paste(f, " + ", "X[[\"", dom, "\"]]*d", l1, 
                sep = "")
        }
    }
    if ((max.level > 1) && (nloc > 1)) {
        for (l1 in 1:(nloc - 1)) {
            for (l2 in (l1 + 1):nloc) {
                aXa <- effNames(c(a, a), c(l1, l2), nloc)
                aXd <- effNames(c(a, d), c(l1, l2), nloc)
                dXa <- effNames(c(d, a), c(l1, l2), nloc)
                dXd <- effNames(c(d, d), c(l1, l2), nloc)
                if (e.unique) {
                  f <- paste(f, " + ", "X[[\"", aXa, "\"]]*a", 
                    l1, "*a", l2, "*ee", sep = "")
                  if (max.dom > 1) {
                    f <- paste(f, " + ", "X[[\"", aXd, "\"]]*a", 
                      l1, "*d", l2, "*ee", sep = "")
                    f <- paste(f, " + ", "X[[\"", dXa, "\"]]*d", 
                      l1, "*a", l2, "*ee", sep = "")
                    f <- paste(f, " + ", "X[[\"", dXd, "\"]]*d", 
                      l1, "*d", l2, "*ee", sep = "")
                  }
                }
                else {
                  f <- paste(f, " + ", "X[[\"", aXa, "\"]]*a", 
                    l1, "*a", l2, "*e", l1, l2, sep = "")
                  if (max.dom > 1) {
                    f <- paste(f, " + ", "X[[\"", aXd, "\"]]*a", 
                      l1, "*d", l2, "*e", l1, l2, sep = "")
                    f <- paste(f, " + ", "X[[\"", dXa, "\"]]*d", 
                      l1, "*a", l2, "*e", l1, l2, sep = "")
                    f <- paste(f, " + ", "X[[\"", dXd, "\"]]*d", 
                      l1, "*d", l2, "*e", l1, l2, sep = "")
                  }
                }
            }
        }
    }
    return(f)
}
