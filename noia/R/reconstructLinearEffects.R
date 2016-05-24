reconstructLinearEffects <-
function (noia.multilinear) 
{
    if (class(noia.multilinear) != "noia.multilinear") {
        stop("Object of class \"multilinear\" expected\n")
    }
    a <- noia::effectsNames[2]
    d <- noia::effectsNames[3]
    e <- noia::effectsNames[4]
    meff <- noia.multilinear$E
    mstd <- noia.multilinear$std.err
    nloc <- noia.multilinear$nloc
    ans.effects <- rep(0, 3^nloc)
    ans.stderr <- rep(0, 3^nloc)
    names(ans.effects) <- effectsNamesGeneral(noia.multilinear$nloc)
    names(ans.stderr) <- effectsNamesGeneral(noia.multilinear$nloc)
    ans.effects[effNames(nloc = nloc)] <- meff[effNames(nloc = nloc)]
    ans.stderr[effNames(nloc = nloc)] <- mstd[effNames(nloc = nloc)]
    for (l1 in 1:nloc) {
        add <- meff[effNames(c(a), c(l1), nloc)]
        dom <- meff[effNames(c(d), c(l1), nloc)]
        std.add <- mstd[effNames(c(a), c(l1), nloc)]
        std.dom <- mstd[effNames(c(d), c(l1), nloc)]
        ans.effects[effNames(c(a), c(l1), nloc)] <- add
        ans.effects[effNames(c(d), c(l1), nloc)] <- dom
        ans.stderr[effNames(c(a), c(l1), nloc)] <- std.add
        ans.stderr[effNames(c(d), c(l1), nloc)] <- std.dom
    }
    if (nloc > 1) {
        for (l1 in 1:(nloc - 1)) {
            for (l2 in (l1 + 1):nloc) {
                a1 <- meff[effNames(c(a), c(l1), nloc)]
                a2 <- meff[effNames(c(a), c(l2), nloc)]
                d1 <- meff[effNames(c(d), c(l1), nloc)]
                d2 <- meff[effNames(c(d), c(l2), nloc)]
                ee <- meff[effNames(c(e, e), c(l1, l2), nloc)]
                cv2.a1 <- ((mstd[effNames(c(a), c(l1), nloc)])/(meff[effNames(c(a), 
                  c(l1), nloc)]))^2
                cv2.a2 <- ((mstd[effNames(c(a), c(l2), nloc)])/(meff[effNames(c(a), 
                  c(l2), nloc)]))^2
                cv2.d1 <- ((mstd[effNames(c(d), c(l1), nloc)])/(meff[effNames(c(d), 
                  c(l1), nloc)]))^2
                cv2.d2 <- ((mstd[effNames(c(d), c(l2), nloc)])/(meff[effNames(c(d), 
                  c(l2), nloc)]))^2
                cv2.ee <- ((mstd[effNames(c(e, e), c(l1, l2), 
                  nloc)])/(meff[effNames(c(e, e), c(l1, l2), 
                  nloc)]))^2
                ans.effects[effNames(c(a, a), c(l1, l2), nloc)] <- a1 * 
                  a2 * ee
                ans.effects[effNames(c(a, d), c(l1, l2), nloc)] <- a1 * 
                  d2 * ee
                ans.effects[effNames(c(d, a), c(l1, l2), nloc)] <- d1 * 
                  a2 * ee
                ans.effects[effNames(c(d, d), c(l1, l2), nloc)] <- d1 * 
                  d2 * ee
                ans.stderr[effNames(c(a, a), c(l1, l2), nloc)] <- sqrt(((a1 * 
                  a2 * ee)^2) * (cv2.a1 + cv2.a2 + cv2.ee + cv2.a1 * 
                  cv2.a2 + cv2.a1 * cv2.ee + cv2.a2 * cv2.ee + 
                  cv2.a1 * cv2.a2 + cv2.ee))
                ans.stderr[effNames(c(a, d), c(l1, l2), nloc)] <- sqrt(((a1 * 
                  d2 * ee)^2) * (cv2.a1 + cv2.d2 + cv2.ee + cv2.a1 * 
                  cv2.d2 + cv2.a1 * cv2.ee + cv2.d2 * cv2.ee + 
                  cv2.a1 * cv2.d2 + cv2.ee))
                ans.stderr[effNames(c(d, a), c(l1, l2), nloc)] <- sqrt(((d1 * 
                  a2 * ee)^2) * (cv2.d1 + cv2.a2 + cv2.ee + cv2.d1 * 
                  cv2.a2 + cv2.d1 * cv2.ee + cv2.a2 * cv2.ee + 
                  cv2.d1 * cv2.a2 + cv2.ee))
                ans.stderr[effNames(c(d, d), c(l1, l2), nloc)] <- sqrt(((d1 * 
                  d2 * ee)^2) * (cv2.d1 + cv2.d2 + cv2.ee + cv2.d1 * 
                  cv2.d2 + cv2.d1 * cv2.ee + cv2.d2 * cv2.ee + 
                  cv2.d1 * cv2.d2 + cv2.ee))
            }
        }
    }
    ans.effects <- ans.effects[colnames(noia.multilinear$smat)]
    ans.stderr <- ans.stderr[colnames(noia.multilinear$smat)]
    return(cbind(ans.effects, ans.stderr))
}
