startingValuesLinear <-
function (noia.linear, max.level = 2, max.dom = 2, e.unique = FALSE, 
    e.init = TRUE) 
{
    if (class(noia.linear) != "noia.linear") {
        stop("Object of class \"noia.linear\" expected\n")
    }
    a <- noia::effectsNames[2]
    d <- noia::effectsNames[3]
    ans <- list()
    nloc <- noia.linear$nloc
    effects <- noia.linear$E
    ans[["R"]] <- effects[effNames(nloc = nloc)]
    for (l1 in 1:nloc) {
        ans[[paste("a", l1, sep = "")]] <- effects[effNames(a, 
            l1, nloc)]
        if (max.dom > 0) {
            ans[[paste("d", l1, sep = "")]] <- effects[effNames(d, 
                l1, nloc)]
        }
    }
    if ((max.level > 1) && (nloc > 1)) {
        if (e.unique) {
            ans[["ee"]] <- 0
        }
        for (l1 in 1:(nloc - 1)) {
            for (l2 in (l1 + 1):nloc) {
                label <- paste("e", l1, l2, sep = "")
                e.unique.factor <- 1
                if (e.unique) {
                  label <- "ee"
                  e.unique.factor <- 1/(nloc * (nloc - 1)/2)
                }
                else {
                  ans[[label]] <- 0
                }
                if (max.dom > 1) {
                  if (e.init) {
                    ans[[label]] <- ans[[label]] + e.unique.factor * 
                      (1/4) * ((effects[effNames(c(a, a), c(l1, 
                      l2), nloc)]/(effects[effNames(a, l1, nloc)] * 
                      effects[effNames(a, l2, nloc)])) + (effects[effNames(c(a, 
                      d), c(l1, l2), nloc)]/(effects[effNames(a, 
                      l1, nloc)] * effects[effNames(d, l2, nloc)])) + 
                      (effects[effNames(c(d, a), c(l1, l2), nloc)]/(effects[effNames(d, 
                        l1, nloc)] * effects[effNames(a, l2, 
                        nloc)])) + (effects[effNames(c(d, d), 
                      c(l1, l2), nloc)]/(effects[effNames(d, 
                      l1, nloc)] * effects[effNames(d, l2, nloc)])))
                  }
                  else {
                    ans[[label]] <- 0
                  }
                }
                else {
                  if (e.init) {
                    ans[[label]] <- ans[[label]] + e.unique.factor * 
                      ((effects[effNames(c(a, a), c(l1, l2), 
                        nloc)]/(effects[effNames(a, l1, nloc)] * 
                        effects[effNames(a, l2, nloc)])))
                  }
                  else {
                    ans[[label]] <- 0
                  }
                }
                if (is.na(ans[[label]])) 
                  ans[[label]] <- 0
            }
        }
    }
    return(ans)
}
